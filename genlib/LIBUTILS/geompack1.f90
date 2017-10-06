!########################################################################################!
!subroutine testpackagegeom1
!subroutine test01
!subroutine test02
!subroutine test03 ( filename )
!subroutine test04
!subroutine test05 ( filename )
!subroutine test06
!subroutine test07
!subroutine test08 ( filename )
!subroutine test09
!subroutine test10
!subroutine test11 ( filename )
!subroutine test12
!subroutine test13 ( filename )
!subroutine test14
!subroutine test15 ( filename )
!subroutine test16
!subroutine test17 ( filename )
!subroutine test18
!subroutine test19
!subroutine pcdec
!subroutine pcds
!subroutine pctri
!subroutine pitri
!subroutine test20
!subroutine test21
!subroutine bedgmv ( nvc, npolg, nvert, maxvc, h, vcl, hvl, pvl, vstart, vnum, &
!subroutine bnsrt2 ( binexp, n, a, map, bin, iwk )
!subroutine cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, &
!subroutine cvdtri ( inter, ldv, nt, vcl, til, tedg, sptr, ierror )
!subroutine delaunay_print ( num_pts, xc, num_tri, nodtri, tnbr )
!subroutine dhpsrt ( k, n, lda, a, map )
!subroutine diam2 ( nvrt, xc, yc, i1, i2, diamsq, ierror )
!subroutine dmat_print ( lda, m, n, a, title )
!subroutine dsftdw ( l, u, k, lda, a, map )
!subroutine dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxho, nvc, npolg, &
!subroutine dsmdf2 ( hflag, nvc, npolg, maxwk, vcl, hvl, pvl, iang, ivrt, &
!subroutine dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv,  &
!subroutine dtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )
!subroutine dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )
!subroutine edght ( a, b, v, n, htsiz, maxedg, hdfree, last, ht, edge, w, ierror )
!subroutine eqdis2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid,  &
!subroutine fndsep ( angac1, xr, yr, nvrt, xc, yc, ivis, theta, nv, iv,  &
!subroutine fndtri ( iedg, mxtr, sflag, tedg, itr, ind, ierror )
!subroutine gtime ( time )
!subroutine holvrt ( nhole, vcl, hvl, pvl, holv )
!subroutine i_swap ( i, j )
!subroutine ihpsrt ( k, n, lda, a, map )
!subroutine imat_print ( lda, m, n, a, title )
!subroutine insed2 ( v, w, npolg, nvert, maxhv, maxpv, vcl, regnum, &
!subroutine insvr2 ( xi, yi, wp, nvc, nvert, maxvc, maxpv, vcl, pvl, &
!subroutine intpg ( nvrt, xc, yc, ctrx, ctry, arpoly, hflag, umdf, wsq, nev, &
!subroutine inttri ( nvrt, xc, yc, h, ibot, costh, sinth, ldv, nvc, ntri,  &
!subroutine isftdw ( l, u, k, lda, a, map )
!subroutine ivec_identity ( n, a )
!subroutine jnhole ( itophv, angspc, angtol, nvc, nvert, maxvc, maxpv,  &
!subroutine lop ( itr, ind, mxedg, top, ldv, vcl, til, tedg, sptr )
!subroutine lufac ( a, lda, n, tol, ipvt, singlr )
!subroutine lusol ( a, lda, n, ipvt, b )
!subroutine mfdec2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid,  &
!subroutine mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )
!subroutine mtredg ( utype, i1, i2, i3, ibndry, nt, til, tedg )
!subroutine prmdf2 ( ipoly, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, listev )
!subroutine ptpolg ( dim, ldv, nv, inc, pgind, vcl, pt, nrml, dtol, inout )
!subroutine randpt ( k, n, seed, axis, nptav, scale, trans, lda, a )
!subroutine resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw,  &
!subroutine rotiar ( n, arr, shift )
!subroutine rotipg ( xeye, yeye, nvrt, xc, yc, ierror )
!subroutine rotpg ( nvrt, xc, yc, i1, i2, ibot, costh, sinth )
!subroutine sepmdf ( angtol, nvrt, xc, yc, arpoly, mean, mdftr, indpvl, &
!subroutine sepshp ( angtol, nvrt, xc, yc, indpvl, iang, i1, i2, wk, ierror )
!subroutine sfdwmf ( l, r, psi, indp, loch )
!subroutine sfupmf ( r, psi, indp, loch )
!subroutine shrnk2 ( nvrt, xc, yc, sdist, nshr, xs, ys, iedge, ierror )
!subroutine spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc,  &
!subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierror )
!subroutine timestamp ( )
!subroutine tmerge ( inter, nbl, ncr, chbl, chcr, ldv, vcl, til, tedg, &
!subroutine trinbr ( nvc, ntri, til, tnbr, htsiz, maxedg, ht, edge, ierror )
!subroutine tripr2 ( nvc, npolg, nvert, maxvc, maxti, maxiw, maxwk, h, vcl,  &
!subroutine trisiz ( ntrid, npolg, hvl, pvl, area, psi, h, indp, loch )
!subroutine trpolg ( nvrt, xc, yc, h, nbc, bndcyc, ldv, nvc, ntri, maxvc,  &
!subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )
!subroutine vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis, ierror )
!subroutine visvrt ( angspc, xeye, yeye, nvis, xc, yc, ivis, maxn, nvsvrt, theta )
!subroutine vornbr ( xeye, yeye, nvrt, xc, yc, nvor, ivor, xvor, yvor, ierror )
!subroutine vpleft ( xc, yc, ivis )
!subroutine vprght ( xc, yc, ivis, ierror )
!subroutine vpscna ( xc, yc, ivis, ierror )
!subroutine vpscnb ( xc, yc, ivis, ierror )
!subroutine vpscnc ( xc, yc, ivis, ierror )
!subroutine vpscnd ( xc, yc, ivis, ierror )
!subroutine walkt2 ( x, y, ntri, vcl, til, tnbr, itri, iedg, ierror )
!subroutine width2 ( nvrt, xc, yc, i1, i2, widsq, ierror )
!subroutine xedge ( mode, xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, xu, yu, intsct )
!subroutine xline ( xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, dv, dw, xu, &
!########################################################################################!



subroutine testpackagegeom1
!***********************************************************************
!
! GEOMPACK_PRB runs the GEOMPACK tests.
!
  implicit none
!
  call timestamp ( ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK_PRB:'
  write ( *, '(a)' ) '  GEOMPACK tests.'

  call test01
  call test02
  call test03 ( 'cmos.in' )
  call test04
  call test05 ( 'annulus.in' )
  call test06
  call test07
  call test08 ( 'annulus.in' )
  call test09
  call test10
  call test11 ( 'annulus.in' )
  call test12
  call test13 ( 'ptpg.in' )
  call test14
  call test15 ( 'shr1.in' )
  call test16
  call test17 ( 'annulus.in' )
  call test18
  call test19
  call test20
  call test21
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  stop
end subroutine

subroutine test01
!
!***********************************************************************
!
!! TEST01 tests ANGLE;
!! TEST01 tests AREAPG;
!! TEST01 tests AREATR.
!
  implicit none
!
  integer, parameter :: n = 9
!
  double precision ang(n)
  double precision angle
  double precision area
  double precision areapg
  double precision areatr
  double precision atr(n)
  integer i
  double precision sum2
  double precision, dimension ( n ) :: xc = (/ &
    5.0D+00, 5.0D+00, 7.0D+00, 9.0D+00, 5.0D+00, &
    1.0D+00, 2.0D+00, 0.0D+00, 3.0D+00 /)
  double precision, dimension ( n ) :: yc = (/ &
    0.0D+00, 2.0D+00, 2.0D+00, 4.0D+00, 8.0D+00, &
    7.0D+00, 5.0D+00, 3.0D+00, 3.0D+00 /)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  AREAPG computes the area of a polygon;'
  write ( *, '(a)') '  AREATR computes the area of a triangle;'
  write ( *, '(a)' ) '  ANGLE computes the polygonal angle at any vertex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of polygonal vertices is N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, XC(I), YC(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,2d15.7)' ) i, xc(i), yc(i)
  end do

  area = areapg ( n, xc, yc )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Area computed directly by AREAPG = ', area

  sum2 = 0.0D+00

  do i = 1, n

   if ( i == 1 ) then
     atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(i+1), yc(i+1) )
     ang(i) = angle ( xc(n), yc(n), xc(i), yc(i), xc(i+1), yc(i+1) )
   else if ( i == n ) then
     atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(1), yc(1) )
     ang(i) = angle ( xc(i-1), yc(i-1), xc(i), yc(i), xc(1), yc(1) )
   else
     atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(i+1), yc(i+1) )
     ang(i) = angle ( xc(i-1), yc(i-1), xc(i), yc(i), xc(i+1), yc(i+1) )
   end if

   sum2 = sum2 + atr(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I, AREATR(I), ANGLE(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,2d15.7)' ) i, atr(i), ang(i)
  end do

  write ( *, '(a,g14.6)' ) &
    'Area computed indirectly by summing AREATR =  ', sum2

  return
end subroutine
subroutine test02
!
!***********************************************************************
!
!! TEST02 tests CMCIRC.
!
  implicit none
!
  integer cmcirc
  integer in
  double precision, parameter :: x0 = 3.0D+00
  double precision, parameter :: y0 =  3.0D+00
  double precision, parameter :: x1 =  5.0D+00
  double precision, parameter :: y1 =  0.0D+00
  double precision, parameter :: x2 =  0.0D+00
  double precision, parameter :: y2 =  5.0D+00
  double precision, parameter :: x3 = -5.0D+00
  double precision, parameter :: y3 =  0.0D+00
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CMCIRC determines if a point lies in, on'
  write ( *, '(a)' ) '  or outside a circle given by 3 points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points defining the circle:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '    X1,Y1 = ', x1, y1
  write ( *, '(a,2g14.6)' ) '    X2,Y2 = ', x2, y2
  write ( *, '(a,2g14.6)' ) '    X3,Y3 = ', x3, y3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point to be tested:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '    X0,Y0 = ', x0, y0

  in = cmcirc ( x0, y0, x1, y1, x2, y2, x3, y3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test results:'

  if ( in == 2 ) then
    write ( *, '(a)' ) '    The three points are collinear.'
  else if ( in == 1 ) then
    write ( *, '(a)' ) '    The point is inside the circle.'
  else if ( in == 0 ) then
    write ( *, '(a)' ) '    The point is on the circle.'
  else if ( in == -1 ) then
    write ( *, '(a)' ) '    The point is outside the circle.'
  end if

  return
end subroutine
subroutine test03 ( filename )
!
!***********************************************************************
!
!! TEST03 tests CVDEC2;
!! TEST03 tests FNDSEP;
!! TEST03 tests INSED2;
!! TEST03 tests INSVR2;
!! TEST03 tests JNHOLE;
!! TEST03 tests MINANG;
!! TEST03 tests RESVRT;
!! TEST03 tests SPDEC2.
!
  implicit none
!
  integer, parameter :: incr = 10000
  integer, parameter :: inunit = 1
  integer, parameter :: maxed = 101
  integer, parameter :: maxhv = 200
  integer, parameter :: maxiw = 900
  integer, parameter :: maxnc = 30
  integer, parameter :: maxnv = 400
  integer, parameter :: maxpv = 1000
  integer, parameter :: maxvc = 500
  integer, parameter :: maxwk = 1500
!
  double precision angmin
  double precision angspc
  double precision angtol
  integer case
  double precision degrees_to_radians
  integer edge(4,maxed)
  character ( len = * ) filename
  integer ht(0:maxed-1)
  integer htsiz
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer icur(maxnc)
  integer ierror
  integer iprt
  integer ivrt(maxnv)
  integer iwk(maxiw)
  integer j
  integer map(maxnc)
  integer msglvl
  integer ncur
  integer nh
  integer nhola
  integer nhole
  integer npolg
  integer nrfv
  integer nsc
  integer nv
  integer nvbc(maxnc)
  integer nvc
  integer nvcin
  integer nvert
  double precision d_pi
  integer prime
  integer pvl(4,maxpv)
  double precision radians_to_degrees
  integer regnum(maxhv)
  character ( len = 20 ) rgname
  double precision tol
  double precision tolin
  double precision vcl(2,maxvc)
  double precision wk(maxwk)
!
  tol = 100.0D+00 * epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!  CASE = -1 or -2: similar, but decompose into simple polygons only
!     and don't obtain convex polygon decomposition
!
  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin, angspc, angtol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) rgname
  write ( *, 700 ) tol,angspc,angtol

  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )

  read ( inunit, * ) case, nvc, ncur, msglvl

  if ( nvc > maxvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRDEC - Error!'
    write ( *, '(a)' ) '  NVC > MAXVC.'
    return
  end if

  if ( ncur > maxnc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRDEC - Error!'
    write ( *, '(a)' ) '  NCUR > MAXNC.'
    return
  end if

  read ( inunit, * ) ( nvbc(i), i = 1, ncur )
  if ( abs(case) == 2 ) then
    read ( inunit, * ) ( icur(i), i = 1, ncur )
  end if
  read ( inunit, * ) ( vcl(1,i), vcl(2,i), i = 1, nvc )
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( abs ( case ) == 1 ) then

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxiw, nvc, npolg, nvert, &
      nhola, regnum, hvl, pvl, iang, iwk, ierror )

  else if ( abs(case) == 2 ) then

    nv = 0

    do i = 1, ncur
      nv = nv + nvbc(i)
    end do

    if ( nv > maxnv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DRDEC - Error!'
      write ( *, '(a)' ) '  NV > MAXNV.'
      return
    end if

    read ( inunit, * ) (ivrt(i),i=1,nv)
    nsc = 0

    do i = 1, nv
      if ( ivrt(i) < 0) then
        nsc = nsc + 1
      end if
    end do

    nsc = nsc / 2

    if ( nsc > maxed ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DRDEC - Error!'
      write ( *, '(a)' ) '  NSC > MAXED.'
      return
    end if

    htsiz = min ( prime ( nsc / 2 ), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, maxiw, &
      npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, iwk, htsiz, nsc, &
      ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRDEC - Error!'
    write ( *, '(a)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole * 2 + nhola
  write ( *, '(a,i6)' ) 'MSGLVL = ', msglvl
  write ( *, '(a,i6)' ) 'NVC = ', nvc
  write ( *,630) (i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) (regnum(i),i=1,npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *,640) nhola,nh,(iwk(i),i=1,nh)

  nrfv = 0
  do i = 1, nvert
    if ( iang(i) > d_pi ( ) + tol ) then
      nrfv = nrfv + 1
    end if
  end do

  angmin = iang(1)
  do i = 2, nvert
    angmin = min ( angmin, iang(i) )
  end do
  angmin = radians_to_degrees ( angmin )

  write (*,710) nvc,npolg,nvert,nhole,nhola,nrfv,angmin
!
!  Obtain simple (and convex) polygon decompositions.
!
  if ( msglvl == 2 ) then
    write ( *,670)
  end if

  call spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc, maxhv, &
    maxpv, maxiw-nh, maxwk, iwk, vcl, regnum, hvl, pvl, iang, iwk(nh+1), wk ,ierror)

  if ( case > 0 ) then

    call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, &
      maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  end if

  if ( msglvl == 2 ) then
    write ( *,670) 0, 0, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRDEC - Error!'
    write ( *, '(a)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after decomposition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NVC = ', nvc
  write ( *,630) (i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,680) npolg, hvl(1:npolg)
  write ( *,650) regnum(1:npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)

  angmin = minval ( iang(1:nvert) )
  angmin = radians_to_degrees ( angmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DECOMP:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin

  630 format ((1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))
  700 format ('input : tol=',d15.7,'   angspc=',f9.3,'   angtol=',f9.3)
  710 format (1x,'initds: nvc=',i7,'   npolg=',i7,'   nvert=',i7, &
    '   nhole=',i7/9x,'nhola=',i7,'   nrfv=',i7,'   angmin=',f9.3 )

  return
end subroutine


subroutine test04
!
!***********************************************************************
!
!! TEST04 tests DIAEDG.
!
  implicit none
!
  integer diaedg
  integer in
  double precision x0
  double precision x1
  double precision x2
  double precision x3
  double precision y0
  double precision y1
  double precision y2
  double precision y3
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DIAEDG determines which diagonal of a'
  write ( *, '(a)' ) '  quadrilateral is to be preferred, based on'
  write ( *, '(a)' ) '  the circumcircle criterion.'

  x0 =  0.0D+00
  y0 =  0.0D+00
  x1 =  5.0D+00
  y1 =  0.0D+00
  x2 =  6.0D+00
  y2 =  1.0D+00
  x3 =  1.0D+00
  y3 =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points defining the quadrilateral:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  P0: X0,Y0 = ', x0, y0
  write ( *, '(a,2g14.6)' ) '  P1: X1,Y1 = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  P2: X2,Y2 = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  P3: X3,Y3 = ', x3, y3

  in = diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIAEDG results:'
  write ( *, '(a)' ) ' '

  if ( in == 1 ) then
    write ( *, '(a)' ) '  Use diagonal P0--P2.'
  else if ( in == -1 ) then
    write ( *, '(a)' ) '  Use diagonal P1--P3.'
  else if ( in == 0 ) then
    write ( *, '(a)' ) '  All 4 points lie on a circle.'
    write ( *, '(a)' ) '  Either diagonal can be used.'
  end if

  return
end subroutine
subroutine test05 ( filename )
!
!***********************************************************************
!
!! TEST05 tests DSMCPR;
!! TEST05 tests DSPGDC;
!! TEST05 tests EDGHT;
!! TEST05 tests HOLVRT.
!
  implicit none
!
  integer, parameter :: incr = 10000
  integer, parameter :: inunit = 1
  integer, parameter :: maxed = 101
  integer, parameter :: maxho = 50
  integer, parameter :: maxnc = 30
  integer, parameter :: maxnv = 500
!
  integer case
  integer edge(4,maxed)
  character ( len = * ) filename
  integer holv(maxho)
  integer ht(0:maxed-1)
  integer htsiz
  integer hvl(maxnc*2)
  integer i
  double precision iang(maxnv*2)
  integer icur(maxnc)
  integer ierror
  integer ivrt(maxnv)
  integer j
  integer map(maxnc)
  integer ncur
  integer nhola
  integer nhole
  integer npolg
  integer nsc
  integer nv
  integer nvbc(maxnc)
  integer nvc
  integer nvert
  integer prime
  integer pvl(4,maxnv*2)
  integer regnum(maxnc*2)
  character ( len = 20 ) rgname
  double precision tolin
  double precision vcl(2,maxnv)
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin
  read ( inunit, * ) case, nvc, ncur

  if ( nvc > maxnv ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Error!'
    write ( *, '(a)' ) '  NVC > MAXNV.'
    return
  end if

  if ( ncur > maxnc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Error!'
    write ( *, '(a)' ) '  NCUR > MAXNC.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)
  if ( case == 2) then
    read ( inunit, * ) icur(1:ncur)
  end if
  read ( inunit, * ) ( vcl(1,i), vcl(2,i), i = 1, nvc )
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxnc*2, maxnv*2, maxho, nvc, &
      npolg, nvert, nhola, regnum, hvl, pvl, iang, holv, ierror )

  else if ( case == 2 ) then

    nv = 0
    do i = 1, ncur
      nv = nv + nvbc(i)
    end do

    if ( nv > maxnv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05 - Error!'
      write ( *, '(a)' ) '  NV > MAXNV.'
      return
    end if

    read ( inunit, * ) ( ivrt(i), i = 1, nv )

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do

    nsc = nsc / 2

    if ( nsc > maxed ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05 - Error!'
      write ( *, '(a)' ) '  NSC > MAXED.'
      return
    end if

    htsiz = min ( prime(nsc/2), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxnc*2, maxnv*2, &
      maxho, npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, holv, &
      htsiz, nsc, ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Error!'
    write ( *, '(a)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'NVC = ', nvc
  write ( *,630) (i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) (regnum(i),i=1,npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *,640) nhola,nhole*2+nhola,(holv(i),i=1,nhole*2+nhola)

  630 format ((1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))

  return
end subroutine

subroutine test06
!
!***********************************************************************
!
!! TEST06 tests DTRIS2;
!
  implicit none
!
  integer, parameter :: maxnp = 10000
  integer, parameter :: maxst = 100
!
  integer a
  integer b
  double precision binexp
  integer c
  integer d
  integer diaedg
  integer i
  integer ierror
  integer ind(maxnp+3)
  integer j
  integer jp1
  integer jp2
  integer k
  integer msglvl
  integer nlo
  integer npt
  integer ntri
  integer stack(maxst)
  integer til(3,maxnp*2+1)
  integer tnbr(3,maxnp*2+1)
  double precision vcl(2,maxnp+3)
!
  msglvl = 0
  binexp = 0.5D+00

  vcl(1,1) =    1.0000000D+00
  vcl(2,1) =    0.0000000D+00
  vcl(1,2) =    0.9238795D+00
  vcl(2,2) =    0.3826834D+00
  vcl(1,3) =    0.7071068D+00
  vcl(2,3) =    0.7071068D+00
  vcl(1,4) =    0.3826834D+00
  vcl(2,4) =    0.9238795D+00
  vcl(1,5) =    0.0000000D+00
  vcl(2,5) =    1.0000000D+00
  vcl(1,6) =  - 0.3826834D+00
  vcl(2,6) =    0.9238795D+00
  vcl(1,7) =  - 0.7071068D+00
  vcl(2,7) =    0.7071068D+00
  vcl(1,8) =  - 0.9238795D+00
  vcl(2,8) =    0.3826834D+00
  vcl(1,9) =  - 1.0000000D+00
  vcl(2,9) =    0.0000000D+00
  vcl(1,10) = - 0.9238795D+00
  vcl(2,10) = - 0.3826834D+00
  vcl(1,11) = - 0.7071068D+00
  vcl(2,11) = - 0.7071068D+00
  vcl(1,12) = - 0.3826834D+00
  vcl(2,12) = - 0.9238795D+00
  vcl(1,13) =   0.0000000D+00
  vcl(2,13) = - 1.0000000D+00
  vcl(1,14) =   0.3826834D+00
  vcl(2,14) = - 0.9238795D+00
  vcl(1,15) =   0.7071068D+00
  vcl(2,15) = - 0.7071068D+00
  vcl(1,16) =   0.9238795D+00
  vcl(2,16) = - 0.3826834D+00
  vcl(1,17) =   0.7500000D+00
  vcl(2,17) =   0.0000000D+00
  vcl(1,18) =   0.6767767D+00
  vcl(2,18) =   0.1767767D+00
  vcl(1,19) =   0.5000000D+00
  vcl(2,19) =   0.2500000D+00
  vcl(1,20) =   0.3232233D+00
  vcl(2,20) =   0.1767767D+00
  vcl(1,21) =   0.2500000D+00
  vcl(2,21) =   0.0000000D+00
  vcl(1,22) =   0.3232233D+00
  vcl(2,22) = - 0.1767767D+00
  vcl(1,23) =   0.5000000D+00
  vcl(2,23) = - 0.2500000D+00
  vcl(1,24) =   0.6767767D+00
  vcl(2,24) = - 0.1767767D+00

  npt = 24

  do i = 1, npt
    ind(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a,i6)' ) '  MSGLVL = ', msglvl
  write ( *, '(a,g14.6)' ) '  BINEXP = ', binexp

  if ( npt > maxnp ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Error!'
    write ( *, '(a)' ) '  NPT > MAXNP.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points to triangulate is ', npt
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The coordinates of the points are:'
  write ( *, '(a)' ) ' '
  do i = 1, npt
    write ( *, '(i5,2f15.7)' ) i, vcl(1,i), vcl(2,i)
  end do

  call dtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  nlo = 0

  do i = 1, ntri

    do j = 1, 3

      k = tnbr(j,i)

      if ( k > i ) then

        jp1 = j + 1
        if ( jp1 > 3) then
          jp1 = 1
        end if

        jp2 = jp1 + 1

        if ( jp2 > 3) then
          jp2 = 1
        end if

        a = til(j,i)
        b = til(jp1,i)
        c = til(jp2,i)

        if ( til(1,k) == b ) then
          d = til(3,k)
        else if ( til(2,k) == b ) then
          d = til(1,k)
        else
          d = til(2,k)
        end if

        if ( diaedg(vcl(1,c),vcl(2,c),vcl(1,a),vcl(2,a),vcl(1,d), &
          vcl(2,d),vcl(1,b),vcl(2,b)) == 1 ) then
          nlo = nlo + 1
        end if

      end if

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'The number of triangles in the triangulation, NTRI = ', ntri
  write ( *, '(a,i6)' ) 'NLO =  ', nlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I, TIL(1:3,I), TNBR(1:3,I)'
  write ( *, '(a)' ) ' '

  do i = 1, ntri
    write ( *, '(7i8)' ) i, (til(j,i),j=1,3), (tnbr(j,i),j=1,3)
  end do


  return
end subroutine

subroutine test07
!
!***********************************************************************
!
!! TEST07 tests DTRIW2;
!
  implicit none
!
  double precision, parameter :: large = 1000.0D+00
  integer, parameter :: maxnp = 10000
  integer, parameter :: maxst = 100
!
  integer a
  integer alg
  integer b
  double precision binexp
  integer c
  integer d
  integer diaedg
  integer i
  integer ierror
  integer ind(maxnp+3)
  integer j
  integer jp1
  integer jp2
  integer k
  integer msglvl
  integer nlo
  integer npt
  integer ntri
  integer stack(maxst)
  integer til(3,maxnp*2+1)
  integer tnbr(3,maxnp*2+1)
  double precision vcl(2,maxnp+3)
!
!  ALG = 
!    2: DTRIW2; 
!    3: DTRIW2 with bounding triangle; 
!    4: DTRIW2 with call to BNSRT2 first.
!  MSGLVL 
!     0: print arrays; 
!     4: also print edges as they are created and swapped.
!
!  I HAVE NO IDEA HOW TO CHOOSE BINEXP
!
  msglvl = 0
  binexp = 0.5D+00

  vcl(1,1) =    1.0000000D+00
  vcl(2,1) =    0.0000000D+00
  vcl(1,2) =    0.9238795D+00
  vcl(2,2) =    0.3826834D+00
  vcl(1,3) =    0.7071068D+00
  vcl(2,3) =    0.7071068D+00
  vcl(1,4) =    0.3826834D+00
  vcl(2,4) =    0.9238795D+00
  vcl(1,5) =    0.0000000D+00
  vcl(2,5) =    1.0000000D+00
  vcl(1,6) =  - 0.3826834D+00
  vcl(2,6) =    0.9238795D+00
  vcl(1,7) =  - 0.7071068D+00
  vcl(2,7) =    0.7071068D+00
  vcl(1,8) =  - 0.9238795D+00
  vcl(2,8) =    0.3826834D+00
  vcl(1,9) =  - 1.0000000D+00
  vcl(2,9) =    0.0000000D+00
  vcl(1,10) = - 0.9238795D+00
  vcl(2,10) = - 0.3826834D+00
  vcl(1,11) = - 0.7071068D+00
  vcl(2,11) = - 0.7071068D+00
  vcl(1,12) = - 0.3826834D+00
  vcl(2,12) = - 0.9238795D+00
  vcl(1,13) =   0.0000000D+00
  vcl(2,13) = - 1.0000000D+00
  vcl(1,14) =   0.3826834D+00
  vcl(2,14) = - 0.9238795D+00
  vcl(1,15) =   0.7071068D+00
  vcl(2,15) = - 0.7071068D+00
  vcl(1,16) =   0.9238795D+00
  vcl(2,16) = - 0.3826834D+00
  vcl(1,17) =   0.7500000D+00
  vcl(2,17) =   0.0000000D+00
  vcl(1,18) =   0.6767767D+00
  vcl(2,18) =   0.1767767D+00
  vcl(1,19) =   0.5000000D+00
  vcl(2,19) =   0.2500000D+00
  vcl(1,20) =   0.3232233D+00
  vcl(2,20) =   0.1767767D+00
  vcl(1,21) =   0.2500000D+00
  vcl(2,21) =   0.0000000D+00
  vcl(1,22) =   0.3232233D+00
  vcl(2,22) = - 0.1767767D+00
  vcl(1,23) =   0.5000000D+00
  vcl(2,23) = - 0.2500000D+00
  vcl(1,24) =   0.6767767D+00
  vcl(2,24) = - 0.1767767D+00

  npt = 24

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a,i6)' ) '  MSGLVL = ', msglvl
  write ( *, '(a,i6)' ) '  NPT =    ', npt
  write ( *, '(a,g14.6)' ) '  BINEXP = ', binexp

  if ( npt > maxnp ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Error!'
    write ( *, '(a)' ) '  NPT > MAXNP.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points to triangulate is ', npt
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The coordinates of the points are:'
  write ( *, '(a)' ) ' '
  do i = 1, npt
    write ( *, '(i5,2f15.7)' ) i, vcl(1,i), vcl(2,i)
  end do

  do alg = 2, 4

    npt = 24

    write ( *, '(a,i6)' ) 'ALG =    ', alg

    if ( alg /= 3 ) then

      do i = 1, npt
        ind(i) = i
      end do

    else

      vcl(1,npt+1) = -large
      vcl(2,npt+1) = -large
      vcl(1,npt+2) = large
      vcl(2,npt+2) = -large
      vcl(1,npt+3) = 0.0D+00
      vcl(2,npt+3) = large
      ind(1) = npt + 1
      ind(2) = npt + 2
      ind(3) = npt + 3
      do i = 1, npt
        ind(i+3) = i
      end do

      npt = npt + 3

    end if

    if ( alg == 4 ) then
      call bnsrt2 ( binexp, npt, vcl, ind, til, tnbr )
    end if

    call dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07 - Error!'
      write ( *, '(a,i6)' ) '  IERROR = ', ierror
      return
    end if

    nlo = 0

    do i = 1, ntri

      do j = 1, 3

        k = tnbr(j,i)

        if ( k > i ) then

          jp1 = j + 1
          if ( jp1 > 3) then
            jp1 = 1
          end if

          jp2 = jp1 + 1

          if ( jp2 > 3) then
            jp2 = 1
          end if

          a = til(j,i)
          b = til(jp1,i)
          c = til(jp2,i)

          if ( til(1,k) == b ) then
            d = til(3,k)
          else if ( til(2,k) == b ) then
            d = til(1,k)
          else
            d = til(2,k)
          end if

          if ( diaedg(vcl(1,c),vcl(2,c),vcl(1,a),vcl(2,a),vcl(1,d), &
               vcl(2,d),vcl(1,b),vcl(2,b)) == 1 ) then
            nlo = nlo + 1
          end if

        end if

      end do

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  NLO =  ', nlo

    call delaunay_print ( npt, vcl, ntri, til, tnbr )

  end do

  return
end subroutine

subroutine test08 ( filename )
!
!***********************************************************************
!
!! TEST08 tests EQDIS2;
!! TEST08 tests INTPG;
!! TEST08 tests MFDEC2;
!! TEST08 tests MMASEP;
!! TEST08 tests SEPMDF;
!! TEST08 tests SEPSHP;
!! TEST08 tests SFDWMF;
!! TEST08 tests SFUPMF;
!! TEST08 tests TRISIZ.
!
  implicit none
!
  integer, parameter :: incr = 10000
  integer, parameter :: inunit = 1
  integer, parameter :: maxed = 101
  integer, parameter :: maxhv = 350
  integer, parameter :: maxiw = 900
  integer, parameter :: maxnc = 30
  integer, parameter :: maxnv = 400
  integer, parameter :: maxpv = 2000
  integer, parameter :: maxvc = 800
  integer, parameter :: maxwk = 1500
!
  double precision angmin
  double precision angspc
  double precision angtol
  double precision area(maxhv)
  integer case
  double precision degrees_to_radians
  double precision dmin
  integer edge(4,maxed)
  character ( len = * ) filename
  double precision h(maxhv)
  logical hflag
  integer ht(0:maxed-1)
  integer htsiz
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer icur(maxnc)
  integer ierror
  integer ivrt(maxnv)
  integer iwk(maxiw)
  integer j
  double precision kappa
  integer map(maxnc)
  integer msglvl
  integer ncur
  integer nh
  integer nhola
  integer nhole
  integer nmin
  integer npolg
  integer nrfv
  integer nsc
  integer ntri(maxhv)
  integer ntrid
  integer ntrie
  integer nv
  integer nvbc(maxnc)
  integer nvc
  integer nvcin
  integer nvert
  double precision d_pi
  integer prime
  double precision psi(maxhv)
  integer pvl(4,maxpv)
  double precision radians_to_degrees
  integer regnum(maxhv)
  character ( len = 20 ) rgname
  double precision tol
  double precision tolin
  double precision umdf2
  double precision vcl(2,maxvc)
  double precision wk(maxwk)
!
  external umdf2
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin,angspc,angtol,kappa,dmin,nmin,ntrid

  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '  Region: ', rgname
  write ( *, 700) tol,angspc,angtol,kappa,dmin,nmin,ntrid

  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )

  hflag = ( kappa >= 0.0D+00 .and. kappa <= 1.0D+00 )
  read ( inunit, * ) case,nvc,ncur,msglvl

  if ( nvc > maxvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error'
    write ( *, '(a)' ) '  NVC > MAXVC.'
    return
  end if

  if ( ncur > maxnc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error'
    write ( *, '(a)' ) '  NCUR > MAXNC.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)
 
  if ( case == 2 ) then
    read ( inunit, * ) icur(1:ncur)
  end if

  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxiw, nvc, npolg, nvert, &
      nhola, regnum, hvl, pvl, iang, iwk, ierror )

  else if ( case == 2 ) then

    nv = sum ( nvbc(1:ncur) )

    if ( nv > maxnv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08 - Error'
      write ( *, '(a)' ) '  NV > MAXNV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do
    nsc = nsc / 2

    if ( nsc > maxed ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08 - Error'
      write ( *, '(a)' ) '  NSC > MAXED.'
      return
    end if

    htsiz = min ( prime ( nsc / 2 ), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, maxiw, &
      npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, iwk, htsiz, nsc, &
      ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole*2 + nhola
  write ( *,670) msglvl
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) (regnum(i),i=1,npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *,640) nhola,nh,(iwk(i),i=1,nh)

  nrfv = 0
  do i = 1, nvert
    if ( iang(i) > d_pi ( ) + tol ) then
      nrfv = nrfv + 1
    end if
  end do

  angmin = minval ( iang(1:nvert) )

  angmin = radians_to_degrees ( angmin )

  write (*,710) nvc,npolg,nvert,nhole,nhola,nrfv,angmin
!
!  Obtain simple and convex polygon decompositions, and print measurements.
!
  if ( msglvl == 2) then
    write ( *,670)
  end if

  call spdec2(angspc,angtol,nvc,npolg,nvert,nhole,nhola,maxvc,maxhv, &
    maxpv,maxiw-nh,maxwk,iwk,vcl,regnum,hvl,pvl,iang,iwk(nh+1),wk,ierror)

  call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, &
    maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  angmin = minval ( iang(1:nvert) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Before call to EQDIS2:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin
!
!  Obtain further convex polygon decomposition based on mesh
!  distribution function, and triangle sizes for the polygons.
!
  write ( *, '(a)' ) 'DEBUG: Call EQDIS2'

  call eqdis2 ( hflag, umdf2, kappa, angspc, angtol, dmin, nmin, ntrid, nvc, &
    npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, &
    iang, area, psi, h, iwk, wk, ierror )

  if ( msglvl == 2 ) then
    write ( *,670) 0, 0, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after decomposition.
!
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,680) npolg, hvl(1:npolg)
  write ( *,650) regnum(1:npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)

  angmin = minval ( iang(1:nvert) )
  angmin = radians_to_degrees ( angmin )

  ntrie = 0
  do i = 1, npolg
    ntri(i) = int ( ntrid * psi(i) * area(i) + 0.5D+00 )
    ntrie = ntrie + ntri(i)
  end do

  write ( *,690) (i,area(i),psi(i),h(i),ntri(i),i=1,npolg)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  After call to EQDIS2:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,i6)' ) '  NTRIE =  ', ntrie
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin

  630 format (/1x,i7/(1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))
  690 format (/(1x,i7,3e15.7,i7))
  700 format ( 'input : tol=',d15.7,'   angspc=',f9.3, &
    '   angtol=',f9.3/9x,'kappa=',f9.3,'   dmin=',f9.3, &
    '   nmin=',i5,'   ntrid=',i7)
  710 format (1x,'initds: nvc=',i7,'   npolg=',i7,'   nvert=',i7, &
    '   nhole=',i7/9x,'nhola=',i7,'   nrfv=',i7,'   angmin=',f9.3 )

  return
end subroutine

subroutine test09
!
!***********************************************************************
!
!! TEST09 tests LRLINE.
!
  implicit none
!
  double precision dv
  integer lr
  integer lrline
  double precision xu
  double precision xv1
  double precision xv2
  double precision yu
  double precision yv1
  double precision yv2
!
  xu = 0.0D+00
  yu = 0.0D+00

  xv1 =  0.0D+00
  yv1 = -1.0D+00
  xv2 =  1.0D+00
  yv2 =  0.0D+00

  dv = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  LRLINE determines if a point is to the right,'
  write ( *, '(a)' ) '  left, or on a directed line that is a directed'
  write ( *, '(a)' ) '  distance away from a directed line from one'
  write ( *, '(a)' ) '  point to another.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The directed base line goes from'
  write ( *, '(2g14.6)' ) xv1, yv1
  write ( *, '(a)' ) '  to'
  write ( *, '(2g14.6)' ) xv2, yv2
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The directed line distance is ', dv
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point to be located is'
  write ( *, '(2g14.6)' ) xu, yu

  lr = lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

  write ( *, '(a)' ) ' '
  if ( lr == -1 ) then
    write ( *, '(a)' ) '  The point is to the right of the line.'
  else if ( lr == 0 ) then
    write ( *, '(a)' ) '  The point is on the line.'
  else if ( lr == 1 ) then
    write ( *, '(a)' ) '  The point is to the left of the line.'
  end if

  return
end subroutine

subroutine test10
!
!***********************************************************************
!
!! TEST10 tests LUFAC;
!! TEST10 tests LUSOL.
!
  implicit none
!
  integer, parameter :: nmax = 20
!
  double precision a(nmax,nmax)
  double precision b(nmax)
  double precision emax
  double precision esum
  integer i
  integer ipvt(nmax)
  integer j
  integer n
  integer seed
  logical singlr
  double precision t
  double precision tol
  real urand
!
  tol = 100.0D+00 * epsilon ( tol )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  LUFAC factors a linear system;'
  write ( *, '(a)' ) '  LUSOL solves a factored linear system.'

  n = 4
  seed = 1952

  b(1:n) = 0.0D+00

  do j = 1, n
    do i = 1, n
      a(i,j) = dble ( urand(seed) ) * 2.0D+00 - 1.0D+00
      b(i) = b(i) + a(i,j)
    end do
  end do

  if ( n <= 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'a, b'
    do i = 1, n
      write ( *, '(5f15.7)' ) (a(i,j),j=1,n),b(i)
    end do
  end if

  call lufac ( a, nmax, n, tol, ipvt, singlr )

  if ( singlr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix is singular'
    return
  end if

  call lusol ( a, nmax, n, ipvt, b )

  if ( n <= 4 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ipvt, lu'
    ipvt(n) = n

    do i = 1, n
      write ( *, '(i5,4f15.7)' ) ipvt(i), (a(i,j),j=1,n)
    end do

  end if

  emax = 0.0D+00
  esum = 0.0D+00
  do i = 1, n
    t = abs ( b(i) - 1.0D+00 )
    emax = max ( emax, t )
    esum = esum + t
  end do

  write ( *,630) (b(i),i=1,min(4,n))
  write ( *,640) emax,esum

  630 format ('  x = ',4f15.7)
  640 format ('  emax,esum = ',2e15.7)

  return
end subroutine

subroutine test11 ( filename )
!
!***********************************************************************
!
!! TEST11 tests DSMDF2;
!! TEST11 tests MDF2;
!! TEST11 tests PRMDF2.
!
  implicit none
!
  integer, parameter :: incr = 10000
  integer, parameter :: inunit = 1
  integer, parameter :: maxed = 101
  integer, parameter :: maxhv = 200
  integer, parameter :: maxiw = 900
  integer, parameter :: maxnc = 30
  integer, parameter :: maxpv = 1000
  integer, parameter :: maxvc = 500
  integer, parameter :: maxwk = 1500
!
  double precision angspc
  double precision angtol
  double precision area(maxhv)
  integer case
  double precision degrees_to_radians
  integer edge(4,maxed)
  double precision edgval(maxpv)
  character ( len = * ) filename
  integer ht(0:maxed-1)
  integer htsiz
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer icur(maxnc)
  integer ierror
  integer ifv
  integer ivrt(maxpv)
  integer iwk(maxiw)
  integer j
  integer map(maxnc)
  double precision mdf2
  integer ncur
  integer nev(maxhv)
  integer nh
  integer nhola
  integer nhole
  integer npolg
  integer nsc
  integer nv
  integer nvbc(maxnc)
  integer nvc
  integer nvert
  integer prime
  integer pvl(4,maxpv)
  integer regnum(maxhv)
  character ( len = 20 ) rgname
  double precision tol
  double precision tolin
  double precision val(maxhv)
  double precision vcl(2,maxvc)
  double precision vrtval(maxvc)
  double precision widsq(maxhv)
  double precision wk(maxwk)
  double precision x
  integer xivrt(maxhv+1)
  double precision y
!
  tol = 100.0D+00 * epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin,angspc,angtol

  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )

  read ( inunit, * ) case,nvc,ncur

  if ( nvc > maxvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error'
    write ( *, '(a)' ) '  NVC > MAXVC.'
    return
  end if

  if ( ncur > maxnc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error'
    write ( *, '(a)' ) '  NCUR > MAXNC.'
    return
  end if

  read ( inunit, * ) (nvbc(i),i=1,ncur)
  if ( case == 2 ) then
    read ( inunit, * ) (icur(i),i=1,ncur)
  end if
  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr(nhole,nvbc,vcl,maxhv,maxpv,maxiw,nvc,npolg,nvert, &
      nhola,regnum,hvl,pvl,iang,iwk, ierror )

  else if ( case == 2 ) then

    nv = 0
    do i = 1, ncur
      nv = nv + nvbc(i)
    end do

    if ( nv > maxpv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11 - Error'
      write ( *, '(a)' ) '  NV > MAXPC.'
      return
    end if

    read ( inunit, * ) (ivrt(i),i=1,nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do
    nsc = nsc/2

    if ( nsc > maxed ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11 - Error'
      write ( *, '(a)' ) '  NSC > MAXED.'
      return
    end if

    htsiz = min(prime(nsc/2),maxed)

    call dspgdc(nvc,vcl,incr,ncur,nvbc,icur,ivrt,maxhv,maxpv,maxiw, &
      npolg,nvert,nhole,nhola,regnum,hvl,pvl,iang,iwk,htsiz,nsc, &
      ht, edge, map, ierror )

  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Obtain simple and convex polygon decompositions.
!
  nh = nhole*2 + nhola

  call spdec2(angspc,angtol,nvc,npolg,nvert,nhole,nhola,maxvc,maxhv, &
    maxpv,maxiw-nh,maxwk,iwk,vcl,regnum,hvl,pvl,iang,iwk(nh+1),wk,ierror)

  call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, &
    maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Initialize data structure for heuristic mesh distribution function
!  and evaluate function at centroid of each convex polygon.
!
  call dsmdf2 ( .true., nvc, npolg, maxwk, vcl, hvl, pvl, iang, &
    ivrt, xivrt, widsq, edgval, vrtval, area, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  do i = 1, npolg

    call prmdf2(i,widsq(i),ivrt,xivrt,edgval,vrtval,nev(i),ifv,iwk)

    if ( nev(i) == 0 ) then

      val(i) = widsq(i)

    else

      nv = xivrt(i+1) - xivrt(i)
      x = 0.0D+00
      y = 0.0D+00
      do j = xivrt(i), xivrt(i+1)-1
        x = x + vcl(1,ivrt(j))
        y = y + vcl(2,ivrt(j))
      end do

      val(i) = 1.0D+00 / mdf2(x/nv,y/nv,widsq(i),nev(i),ifv,iwk,ivrt, &
        edgval,vrtval,vcl)

    end if

  end do

  close ( unit = inunit )
!
!  Print arrays from calls to 3 mdf routines.
!
  write ( *,630) npolg,(i,xivrt(i),area(i),widsq(i),val(i),nev(i),i=1,npolg)
  write ( *,640) nvert,nvc,(i,ivrt(i),edgval(i),vrtval(i),i=1,nvc)
  write ( *,650) (i,ivrt(i),edgval(i),i=nvc+1,nvert)

  630 format (/1x,i7/(1x,2i7,3f15.7,i7))
  640 format (/1x,2i7/(1x,2i7,2f15.7))
  650 format (1x,2i7,f15.7)

  return
end subroutine


subroutine test12
!
!***********************************************************************
!
!! TEST12 tests PRIME.
!
  implicit none
!
  integer i
  integer prime
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  PRIME finds the smallest prime bigger than'
  write ( *, '(a)' ) '  a given value.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, PRIME(I)'
  write ( *, '(a)' ) ' '

  do i = 100, 500, 100
    write ( *, '(2i10)' ) i, prime ( i )
  end do

  return
end subroutine

subroutine test13 ( filename )
!
!***********************************************************************
!
!! TEST13 tests PTPOLG.
!
  implicit none
!
  integer, parameter :: inunit = 1
  integer, parameter :: maxn = 100
!
  double precision a
  double precision b
  double precision c
  integer dim
  double precision dtol
  character ( len = * ) filename
  integer i
  integer inout
  integer j
  integer n
  double precision nrml(3)
  integer pgind(0:maxn)
  double precision pt(3)
  double precision tol
  double precision vcl(3,maxn)
!
  tol = 100.0D+00 * epsilon ( tol )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  dtol = 10.0D+00 * tol
  read ( inunit, * ) n,dim,a,b
  if ( n < 3 ) then
    return
  end if

  if ( n > maxn ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRPTPG - Error'
    write ( *, '(a)' ) '  N > MAXN.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'DIM = ', dim
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) ' '

  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,n)

  pgind(0) = n
  do i = 1, n
    pgind(i) = i
    if ( dim == 3 ) then
      vcl(3,i) = a * vcl(1,i) + b * vcl(2,i)
    end if
    write ( *, '(i5,3f15.7)' ) i,(vcl(j,i),j=1,dim)
  end do

  if ( dim == 3 ) then
    c = sqrt ( 1.0D+00 + a**2 + b**2)
    nrml(1) = - a / c
    nrml(2) = - b / c
    nrml(3) = 1.0D+00 / c
  end if

  write ( *, '(a)' ) ' '

   20 continue

  read ( inunit, *, end=30 ) pt(1),pt(2)

  if ( dim == 3 ) then
    pt(3) = a * pt(1) + b * pt(2)
  end if

  call ptpolg ( dim, 3, n, 1, pgind, vcl, pt, nrml, dtol, inout )

  write ( *,630) inout,(pt(i),i=1,dim)
  go to 20

   30 continue

  close ( unit = inunit )

  630 format (1x,'inout=',i3,3x,'pt=',3f15.7)

  return
end subroutine

subroutine test14
!
!***********************************************************************
!
!! TEST14 tests ROTIAR.
!
  implicit none
!
  integer, parameter :: maxn = 50
!
  integer a(maxn)
  integer i
  integer n
  integer shift
!
  n = 10
  shift = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  ROTIAR "rotates" an array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Using N = ', n
  write ( *, '(a,i6)' ) '  SHIFT = ', shift

  do i = 1, n
    a(i) = i
  end do

  write ( *, '(10i6)' ) a(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shifted array:'
  write ( *, '(a)' ) ' '

  call rotiar ( n, a, shift )

  write ( *, '(10i6)' ) a(1:n)

  return
end
subroutine test15 ( filename )
!
!***********************************************************************
!
!! TEST15 tests DIAM2;
!! TEST15 tests SHRNK2;
!! TEST15 tests WIDTH2.
!
  implicit none
!
  integer, parameter :: inunit = 1
  integer, parameter :: maxn = 100
!
  double precision diamsq
  character ( len = * ) filename
  integer i
  integer i1
  integer i2
  integer iedge(0:maxn)
  integer ierror
  integer n
  integer nshr
  double precision sdist(0:maxn-1)
  double precision widsq
  double precision xc(0:maxn)
  double precision xs(0:maxn)
  double precision yc(0:maxn)
  double precision ys(0:maxn)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  Reading data file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, * ) n

  if ( n > maxn ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  Maximum legal N is ', maxn
    return
  end if

  read ( inunit, * ) (xc(i),yc(i),sdist(i),i=0,n-1)

  close ( unit = inunit )

  xc(n) = xc(0)
  yc(n) = yc(0)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(2f15.7)' ) (xc(i),yc(i),i=0,n)
  xs(0) = 0.0D+00
  ys(0) = 0.0D+00

  call shrnk2 ( n, xc, yc, sdist, nshr, xs, ys, iedge ,ierror)

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  SHRNK2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NSHR = ', n
  write ( *, '(2f15.7)' ) (xs(i),ys(i),i=0,nshr)

  call diam2 ( n, xc(1), yc(1), i1, i2, diamsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  DIAM2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,diamsq

  call width2 ( n, xc(1), yc(1), i1, i2, widsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  WIDTH2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,widsq

  if ( nshr < 3 ) then
    return
  end if

  call diam2 ( nshr, xs(1), ys(1), i1, i2, diamsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  DIAM2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,diamsq

  call width2 ( nshr, xs(1), ys(1), i1, i2, widsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  WIDTH2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,widsq

  return
end subroutine

subroutine test16
!
!***********************************************************************
!
!! TEST16 tests DHPSRT;
!! TEST16 tests IHPSRT;
!! TEST16 tests RANDPT.
!
  implicit none
!
  integer, parameter :: maxk = 4
  integer, parameter :: maxn = 100
!
  integer axis
  double precision da(maxk,maxn)
  integer i
  integer ia(maxk,maxn)
  logical iflag
  integer j
  integer k
  integer map(maxn)
  integer n
  integer nptav
  double precision, dimension ( maxk ) :: scale = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  integer seed
  double precision, dimension ( maxk ) :: trans = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'

  k = 2
  n = 10
  seed = 1952
  axis = 0
  nptav = 0

  iflag = (k < 0)
  k = abs(k)

  if ( k < 1 .or. k > maxk ) then
    return
  end if

  if ( n < 1 .or. n > maxn ) then
    return
  end if

  do j = 1, n
    map(j) = j
  end do

  call randpt ( k, n, seed, axis, nptav, scale, trans, maxk, da )

  if ( iflag ) then

    do j = 1, n
      do i = 1, k
        ia(i,j) = int ( n * da(i,j) )
      end do
    end do

    do j = 1, n
      write ( *, '(5i5)' ) j,(ia(i,j),i=1,k)
    end do

    call ihpsrt ( k, n, maxk, ia, map )

    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(5i5)' ) map(j),(ia(i,map(j)),i=1,k)
    end do

  else

    do j = 1, n
      write ( *, '(i5,4f15.7)' ) j,(da(i,j),i=1,k)
    end do

    call dhpsrt ( k, n, maxk, da, map )

    write ( *,'(a)' ) ' '
    do j = 1, n
     write ( *, '(i5,4f15.7)' ) map(j), (da(i,map(j)),i=1,k)
    end do

  end if

  return
end subroutine

subroutine test17 ( filename )
!
!***********************************************************************
!
!! TEST17 tests BEDGMV;
!! TEST17 tests CVDTRI;
!! TEST17 tests FNDTRI;
!! TEST17 tests INTTRI;
!! TEST17 tests LOP;
!! TEST17 tests MTREDG;
!! TEST17 tests ROTPG;
!! TEST17 tests TMERGE;
!! TEST17 tests TRINBR;
!! TEST17 tests TRIPR2;
!! TEST17 tests TRPOLG.
!
  implicit none
!
  integer, parameter :: incr = 10000
  integer, parameter :: inunit = 1
  integer, parameter :: maxed = 1201
  integer, parameter :: maxhv = 350
  integer, parameter :: maxiw = 900
  integer, parameter :: maxnc = 30
  integer, parameter :: maxnv = 400
  integer, parameter :: maxpv = 2000
  integer, parameter :: maxti = 8000
  integer, parameter :: maxvc = 5000
  integer, parameter :: maxwk = 1500
  integer, parameter :: nfreq = 12
!
  integer a
  integer afreq(0:nfreq-1)
  double precision anga
  double precision angb
  double precision angc
  double precision angle
  double precision angmax
  double precision angmin
  double precision angspc
  double precision angtol
  double precision area(maxhv)
  integer b
  integer c
  integer case
  double precision degrees_to_radians
  double precision delta
  double precision dmin
  integer edge(4,maxed)
  character ( len = * ) filename
  double precision h(maxhv)
  logical hflag
  integer ht(0:maxed-1)
  integer htsiz
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer icur(maxnc)
  integer ierror
  integer ivrt(maxnv)
  integer iwk(maxiw)
  integer j
  double precision kappa
  integer map(maxnc)
  integer msglvl
  integer ncur
  integer nh
  integer nhola
  integer nhole
  integer nmin
  integer npolg
  integer nrfv
  integer nsc
  integer ntri
  integer ntrid
  integer nv
  integer nvbc(maxnc)
  integer nvc
  integer nvcin
  integer nvert
  double precision d_pi
  integer prime
  double precision psi(maxhv)
  integer pvl(4,maxpv)
  double precision radians_to_degrees
  integer regnum(maxhv)
  character ( len = 20 ) rgname
  integer til(3,maxti)
  integer tnbr(3,maxti)
  double precision tol
  double precision tolin
  integer tstart(maxhv)
  double precision, external :: umdf2
  double precision vcl(2,maxvc)
  integer vnum(maxpv)
  integer vstart(maxpv)
  double precision wk(maxwk)
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin, angspc, angtol, kappa, dmin, nmin, ntrid

  write ( *, '(a)' ) ' '
  write ( *, 710 ) rgname, tol, angspc, angtol, kappa, dmin, nmin, ntrid
  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )
  hflag = ( kappa >= 0.0D+00 .and. kappa <= 1.0D+00 )

  read ( inunit, * ) case, nvc, ncur, msglvl

  if ( nvc > maxvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error'
    write ( *, '(a)' ) '  NVC > MAXVC.'
    return
  end if

  if ( ncur > maxnc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error'
    write ( *, '(a)' ) '  NCUR > MAXNC.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)

  if ( case == 2 ) then
    read ( inunit, * ) icur(1:ncur)
  end if

  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call DSMCPR or DSPGDC to set the data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxiw, nvc, npolg, nvert, &
      nhola, regnum, hvl, pvl, iang, iwk, ierror )

  else if ( case == 2 ) then

    nv = sum ( nvbc(1:ncur) )

    if ( nv > maxnv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17 - Error'
      write ( *, '(a)' ) '  NV > MAXNV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do

    nsc = nsc / 2

    if ( nsc > maxed ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17 - Error'
      write ( *, '(a)' ) '  NSC > MAXED.'
      return
    end if

    htsiz = min ( prime(nsc/2), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, maxiw, &
      npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, iwk, htsiz, nsc, &
      ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out the data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole * 2 + nhola
  write ( *, '(a,i6)' ) 'MSGLVL = ', msglvl
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) regnum(1:npolg)
  write ( *,660) nvert, (i,(pvl(j,i),j=1,4), iang(i), i = 1,nvert)
  write ( *,640) nhola, nh, iwk(1:nh)

  nrfv = 0
  do i = 1, nvert
    if ( iang(i) > d_pi ( ) + tol ) then
      nrfv = nrfv + 1
    end if
  end do

  angmin = minval ( iang(1:nvert) )
  angmin = radians_to_degrees ( angmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17: INITDS:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,i6)' ) '  NHOLE =  ', nhole
  write ( *, '(a,i6)' ) '  NHOLA =  ', nhola
  write ( *, '(a,i6)' ) '  NRFV =   ', nrfv
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin
!
!  Obtain simple and convex polygon decompositions, and print measurements.
!
  if ( msglvl == 2 ) then
    write ( *,670)
  end if

  call spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc, maxhv, &
    maxpv, maxiw-nh, maxwk, iwk, vcl, regnum, hvl, pvl, iang, iwk(nh+1), wk ,ierror)

  call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, &
    maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  angmin = minval ( iang(1:nvert) )

  angmin = radians_to_degrees ( angmin )

  write (*,730) nvc, npolg, nvert, angmin
!
!  Obtain further convex polygon decomposition based on mesh
!  distribution function, and triangle sizes for the polygons.
!  Then print measurements.
!
  call eqdis2 ( hflag, umdf2, kappa, angspc, angtol, dmin, nmin, ntrid, nvc, &
    npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, &
    iang, area, psi, h, iwk, wk, ierror )

  if ( msglvl == 2 ) then
    write ( *,670) 0, 0, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,680) npolg, hvl(1:npolg)
  write ( *,650) regnum(1:npolg)

  angmin = radians_to_degrees ( minval ( iang(1:nvert) ) )

  write (*,740) nvc, npolg, nvert, angmin
  htsiz = min ( prime(nvcin*2), maxed )
  nvcin = nvc
!
!  Triangulate each convex polygon in the decomposition, according to
!  the mesh spacings in the H array.
!
  call tripr2 ( nvc, npolg, nvert, maxvc, maxti, maxiw, maxwk, h, vcl, hvl, &
    pvl, iang, ntri, til, vstart, vnum, tstart, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Compute TNBR array. Then print arrays and measurements.
!
  call trinbr ( nvc, ntri, til, tnbr, htsiz, maxed, ht, edge ,ierror)

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  write ( *,690) nvert,(i,(pvl(j,i),j=1,4),iang(i),vstart(i),vnum(i),i=1,nvert)
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,650) tstart(1:npolg)
  write ( *,700) ntri,(i,(til(j,i),j=1,3),(tnbr(j,i),j=1,3),i=1,ntri)

  angmin = d_pi ( )
  angmax = 0.0D+00
  delta = d_pi ( ) / dble(nfreq)
  afreq(0:nfreq-1) = 0

  do i = 1, ntri

    a = til(1,i)
    b = til(2,i)
    c = til(3,i)

    anga = angle(vcl(1,c),vcl(2,c),vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b))

    angb = angle(vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),vcl(1,c),vcl(2,c))

    angc = d_pi ( ) - anga - angb
    angmin = min ( angmin, anga, angb, angc )
    angmax = max ( angmax, anga, angb, angc )

    if ( abs ( anga - 0.5D+00 * d_pi ( ) ) <= tol ) then
      anga = 0.5D+00 * d_pi ( ) - 2.0D+00 * tol
    end if

    j = min ( nfreq - 1, int ( ( anga + tol ) / delta ) )
    afreq(j) = afreq(j) + 1

    if ( abs ( angb - 0.5D+00 * d_pi ( ) ) <= tol ) then
      angb = 0.5D+00 * d_pi ( ) - 2.0D+00 * tol
    end if

    j = min ( nfreq - 1, int ( ( angb + tol ) / delta ) )
    afreq(j) = afreq(j) + 1

    if ( abs ( angc - 0.5D+00 * d_pi ( ) ) <= tol ) then
      angc = 0.5D+00 * d_pi ( ) - 2.0D+00 * tol
    end if

    j = min ( nfreq - 1, int ( ( angc + tol ) / delta ) )
!
!  Having some odd cases where J comes out NEGATIVE.
!
    if ( j < 0 ) then
      j = 0
    end if

    afreq(j) = afreq(j) + 1

  end do

  angmin = radians_to_degrees ( angmin )
  angmax = radians_to_degrees ( angmax )
  write (*,750) nvc,ntri,angmin,angmax
  write (*,760) (i*180.0D+00/dble(nfreq),afreq(i)/dble(3*ntri),i=0,nfreq-1)

  630 format (/1x,i7/(1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))
  690 format (/1x,i7/(1x,5i7,f15.7,2i7))
  700 format (/1x,i7/(1x,7i7))
  710 format (1x,a20/1x,'input : tol=',d15.7,'   angspc=',f9.3, &
        '   angtol=',f9.3/9x,'kappa=',f9.3,'   dmin=',f9.3, &
        '   nmin=',i5,'   ntrid=',i7)
  730 format (1x,'decomp: nvc=',i7,'   npolg=',i7,'   nvert=',i7/ &
        9x,'angmin=',f9.3 )
  740 format (1x,'eqdist: nvc=',i7,'   npolg=',i7,'   nvert=',i7/ &
        9x,'angmin=',f9.3 )
  750 format (1x,'triang: nvc=',i7,'   ntri=',i7 / &
        9x,'angmin=',f9.3,'   angmax=',f9.3 )
  760 format (1x,'relative frequency of triangle angles'/4(f8.1,f10.5))

  return
end subroutine


subroutine test18
!
!***********************************************************************
!
!! TEST18 tests VISVRT;
!! TEST18 tests VORNBR.
!
  implicit none
!
  integer, parameter :: maxn = 200
!
  double precision angle
  double precision angspc
  double precision degrees_to_radians
  integer i
  integer ierror
  integer ivert
  integer ivis(0:maxn)
  integer ivor(0:maxn)
  integer j
  integer maxnv
  integer n
  integer nvis
  integer nvor
  integer nvrt
  integer nvsvrt
  double precision phi
  double precision theta(0:maxn)
  double precision x(maxn)
  double precision xc(0:maxn)
  double precision xeye
  double precision xvor(0:maxn)
  double precision y(maxn)
  double precision yc(0:maxn)
  double precision yeye
  double precision yvor(0:maxn)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  Test being SKIPPED FOR NOW.'
  return

  read ( *, * ) n

  if ( n < 3 ) then
    return
  end if

  if ( n > maxn ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST18 - Error'
    write ( *, '(a)' ) '  N > MAXN.'
    return
  end if

  read ( *, * ) ( x(i), y(i), i = 1, n )
  read ( *, * ) ivert
  read ( *, * ) angspc
!
!  1 <= IVERT <= N is index of polygon vertex for eyepoint.
!  ANGSPC is angle spacing parameter in degrees.
!
  angspc = degrees_to_radians ( angspc )
  xeye = x(ivert)
  yeye = y(ivert)
  nvrt = n - 2

  j = -1
  do i = ivert+1, n
    j = j + 1
    xc(j) = x(i)
    yc(j) = y(i)
  end do

  do i = 1, ivert-1
    j = j + 1
    xc(j) = x(i)
    yc(j) = y(i)
  end do

  write ( *, '(a,2g14.6)' ) '(XEYE,YEYE)= ', xeye, yeye
  write ( *,630) nvrt,(xc(i),yc(i),i=0,nvrt)

  call vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis,ierror )

  write ( *,640) nvis,(xc(i),yc(i),ivis(i),i=0,nvis)

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST18 - Fatal error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  phi = angle ( xc(nvis), yc(nvis), xeye, yeye, xc(0), yc(0) )
  maxnv = nvis + int ( phi / angspc )

  if ( maxnv > maxn ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST18 - Error'
    write ( *, '(a)' ) '  MAXNV > MAXN.'
    return
  end if

  call visvrt ( angspc, xeye, yeye, nvis, xc, yc, ivis, maxnv, nvsvrt, theta )

  write ( *,650) nvsvrt,(xc(i),yc(i),ivis(i),theta(i),i=0,nvsvrt)
  nvor = -1

  call vornbr ( xeye, yeye, nvsvrt, xc, yc, nvor, ivor, xvor, yvor, ierror )

  write ( *,640) nvor,(xvor(i),yvor(i),ivor(i),i=0,nvor)

  630 format (/1x,i5/(1x,2f15.7))
  640 format (/1x,i5/(1x,2f15.7,i5))
  650 format (/1x,i5/(1x,2f15.7,i5,f15.7))

  return
end subroutine


subroutine test19
!
!***********************************************************************
!
!! TEST19 tests ROTIPG;
!! TEST19 tests VISPOL.
!
  implicit none
!
  integer, parameter :: maxn = 200
!
  integer i
  integer ierror
  integer ivert
  integer ivis(0:maxn+1)
  integer j
  integer n
  integer nvis
  integer nvrt
  integer vptype
  double precision x(maxn)
  double precision xc(0:maxn+1)
  double precision xeye
  double precision y(maxn)
  double precision yc(0:maxn+1)
  double precision yeye
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  TEST BEING SKIPPED FOR NOW.'
  return

  read ( *, * ) n
  if ( n < 3 ) then
    return
  end if

  if ( n > maxn ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST19 - Error'
    write ( *, '(a)' ) '  N > MAXN.'
    return
  end if

  read ( *,*) (x(i),y(i),i=1,n)
  read ( *,*) ivert,vptype,xeye,yeye
!
!  1 <= IVERT <= N is index of polygon vertex or IVERT <= 0 if
!  ROTIPG is to be called. BVTYPE = 0, 1, or 2 for boundary,
!  interior, or blocked exterior viewpoint. (XEYE,YEYE) is needed
!  for non-boundary viewpoints only, and must be visible from
!  (X(IVERT),Y(IVERT)) if IVERT > 0.
!
  if ( vptype == 0 ) then

    xeye = x(ivert)
    yeye = y(ivert)
    nvrt = n - 2

    j = -1
    do i = ivert+1, n
      j = j + 1
      xc(j) = x(i)
      yc(j) = y(i)
    end do

    do i = 1, ivert-1
      j = j + 1
      xc(j) = x(i)
      yc(j) = y(i)
    end do

  else if ( vptype == 1 ) then

    nvrt = n

    if ( ivert <= 0 ) then

      do i = 1, n
        xc(i) = x(i)
        yc(i) = y(i)
      end do

      xc(0) = xc(n)
      yc(0) = yc(n)

      call rotipg ( xeye, yeye, nvrt, xc, yc ,ierror)

    else

      j = -1

      do i = ivert, n
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

      do i = 1, ivert
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

    end if

  else

    nvrt = n

    if ( ivert <= 0 ) then

      do i = 1, n
        xc(n-i) = x(i)
        yc(n-i) = y(i)
      end do

      xc(n) = xc(0)
      yc(n) = yc(0)

      call rotipg ( xeye, yeye, nvrt, xc, yc ,ierror)

    else

      j = -1

      do i = ivert, 1, -1
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

      do i = n, ivert, -1
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

    end if

  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST19 - Error.'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '(XEYE,YEYE)=', xeye,yeye
  write ( *,630) nvrt,(xc(i),yc(i),i=0,nvrt)

  call vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis ,ierror)

  write ( *,640) nvis,(xc(i),yc(i),ivis(i),i=0,nvis)

  630 format (/1x,i5/(1x,2f15.7))
  640 format (/1x,i5/(1x,2f15.7,i5))

  return
end subroutine


subroutine pcdec
!
!***********************************************************************
!
!! PCDEC generates PICTEX commands from DRDEC or DREQD output.
!
  implicit none
!
  integer, parameter :: edgv = 4
  integer, parameter :: loc = 1
  integer, parameter :: maxho = 50
  integer, parameter :: maxnc = 30
  integer, parameter :: maxnv = 500
  integer, parameter :: succ = 3
!
  integer holv(maxho)
  integer hvl(maxnc*2)
  integer i
  double precision iang(maxnv*2)
  integer j
  integer j1
  integer l
  integer l1
  logical mid
  integer msglvl
  integer nh
  integer nhola
  integer nhole
  integer npolg
  integer nvc
  integer nvert
  integer pvl(4,maxnv*2)
  integer regnum(maxnc*2)
  integer t
  double precision vcl(2,maxnv)
  double precision x1
  double precision x2
  double precision xmax
  double precision xmin
  double precision y1
  double precision y2
  double precision ymax
  double precision ymin
!
  read ( *,*) msglvl
  if ( msglvl /= 2) then
    return
  end if
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *,*) npolg,nhole
  read ( *,*) (hvl(i),i=1,npolg+nhole)
  read ( *,*) (regnum(i),i=1,npolg)
  read ( *,*) nvert
  read ( *,*) (t,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  read ( *,*) nhola,nh
  read ( *,*) (holv(i),i=1,nh)
 
  xmin = vcl(1,1)
  ymin = vcl(2,1)
  xmax = xmin
  ymax = ymin
  do i = 2, nvc
    xmin = min(xmin,vcl(1,i))
    xmax = max(xmax,vcl(1,i))
    ymin = min(ymin,vcl(2,i))
    ymax = max(ymax,vcl(2,i))
  end do

  write ( *, '(a)' ) '\\begin{figure}[htbp]'
  write ( *, '(a)' ) '\\begin{center}'
  write ( *, '(a)' ) '\\mbox{\\beginpicture'
  write ( *, '(a)' ) '\\setcoordinatesystem units <0.2cm,0.2cm>'
  write ( *,610) xmin,xmax,ymin,ymax
  write ( *, '(a)' ) '\\setlinear'

  do i = 1, npolg+nhole

    mid = .false.
    j = hvl(i)
    l = pvl(loc,j)

   20   continue

    j1 = pvl(succ,j)
    l1 = pvl(loc,j1)

    if ( pvl(edgv,j) < j ) then

      if ( .not. mid ) then
        write ( *, '(a)' ) '\\plot'
        mid = .true.
        write ( *, '(2f15.7)' ) vcl(1,l),vcl(2,l)
      end if

      write ( *, '(2f15.7)' ) vcl(1,l1),vcl(2,l1)

    else if ( mid ) then

      write ( *, '(a)' ) '/'
      mid = .false.

    end if

    j = j1
    l = l1

    if ( j /= hvl(i) ) then
      go to 20
    end if

    if ( mid ) then
      write ( *, '(a)' ) '/'
    end if

  end do

   40 continue

   read ( *,*) i,j,x1,y1,x2,y2
   if ( i == 0 ) then
     go to 50
   end if
   write ( *,630) x1,y1,x2,y2
   go to 40

   50 continue

  write ( *, '(a)' ) '\\endpicture}'
  write ( *, '(a)' ) '\\end{center}'
  write ( *, '(a)' ) '\\end{figure}'

  610 format ('\\setplotarea x from ',f9.4,' to ',f9.4,', y from ',f9.4, &
    ' to ',f9.4)
  630 format ('\\plot ',4f15.7,' /')

  return
end subroutine



subroutine pcds
!
!***********************************************************************
!
!! PCDS generates PICTEX commands from DRDS output.
!
  implicit none
!
  integer, parameter :: edgv = 4
  integer, parameter :: loc = 1
  integer, parameter :: maxho = 50
  integer, parameter :: maxnc = 30
  integer, parameter :: maxnv = 500
  integer, parameter :: succ = 3
!
  integer holv(maxho)
  integer hvl(maxnc*2)
  integer i
  double precision iang(maxnv*2)
  integer j
  integer j1
  integer l
  integer l1
  logical mid
  integer nh
  integer nhola
  integer nhole
  integer npolg
  integer nvc
  integer nvert
  integer pvl(4,maxnv*2)
  integer regnum(maxnc*2)
  integer t
  double precision vcl(2,maxnv)
  double precision xmax
  double precision xmin
  double precision ymax
  double precision ymin
!
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *,*) npolg,nhole
  read ( *,*) (hvl(i),i=1,npolg+nhole)
  read ( *,*) (regnum(i),i=1,npolg)
  read ( *,*) nvert
  read ( *,*) (t,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  read ( *,*) nhola,nh
  read ( *,*) (holv(i),i=1,nh)

  xmin = vcl(1,1)
  ymin = vcl(2,1)
  xmax = xmin
  ymax = ymin
  do i = 2, nvc
    xmin = min(xmin,vcl(1,i))
    xmax = max(xmax,vcl(1,i))
    ymin = min(ymin,vcl(2,i))
    ymax = max(ymax,vcl(2,i))
  end do

  write ( *, '(a)' ) '\\begin{figure}[htbp]'
  write ( *, '(a)' ) '\\begin{center}'
  write ( *, '(a)' ) '\\mbox{\\beginpicture'
  write ( *, '(a)' ) '\\setcoordinatesystem units <0.2cm,0.2cm>'
  write ( *, 610) xmin,xmax,ymin,ymax
  write ( *, '(a)' ) '\\setlinear'

  do i = 1, npolg+nhole

    mid = .false.
    j = hvl(i)
    l = pvl(loc,j)

   20   continue

    j1 = pvl(succ,j)
    l1 = pvl(loc,j1)

    if ( pvl(edgv,j) < j ) then

      if ( .not. mid ) then
        write ( *, '(a)' ) '\\plot'
        mid = .true.
        write ( *, '(2f15.7)' ) vcl(1,l),vcl(2,l)
      end if

      write ( *, '(2f15.7)' ) vcl(1,l1),vcl(2,l1)

    else if ( mid ) then

      write ( *, '(a)' ) '/'
      mid = .false.

    end if

    j = j1
    l = l1

    if ( j /= hvl(i) ) then
      go to 20
    end if

    if ( mid ) then
      write ( *, '(a)' ) '/'
    end if

  end do

  write ( *, '(a)' ) '\\endpicture}'
  write ( *, '(a)' ) '\\end{center}'
  write ( *, '(a)' ) '\\end{figure}'

  610 format ('\\setplotarea x from ',f9.4,' to ',f9.4,', y from ',f9.4, &
    ' to ',f9.4)

  return
end subroutine



subroutine pctri
!
!***********************************************************************
!
!! PCTRI generates PICTEX commands from DRTRI output.
!
  implicit none
!
  integer, parameter :: maxti = 8000
  integer, parameter :: maxvc = 5000
!
  integer a
  integer b
  integer i
  integer j
  integer jp1
  integer msglvl
  integer nh
  integer nhola
  integer nhole
  integer npolg
  integer ntri
  integer nvc
  integer nvcin
  integer nvert
  integer t
  integer til(3,maxti)
  integer tnbr(3,maxti)
  double precision vcl(2,maxvc)
  double precision x
  double precision xmax
  double precision xmin
  double precision ymax
  double precision ymin
!
  read ( *,*,end=40) xmin,xmax,ymin,ymax
  read ( *,*) msglvl
  if ( msglvl /= 0) then
    return
  end if
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *,*) npolg,nhole
  read ( *,*) (t,i=1,npolg+nhole)
  read ( *,*) (t,i=1,npolg)
  read ( *,*) nvert
  read ( *,*) ((t,j=1,5),x,i=1,nvert)
  read ( *,*) nhola,nh
  read ( *,*) (t,i=1,nh)
  nvcin = nvc
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *,*) npolg
  read ( *,*) (t,i=1,npolg)
  read ( *,*) (t,i=1,npolg)
  read ( *,*) nvert
  read ( *,*) ((t,j=1,5),x,t,t,i=1,nvert)
  nvcin = nvc
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *,*) (t,i=1,npolg)
  read ( *,*) ntri
  read ( *,*) (t,(til(j,i),j=1,3),(tnbr(j,i),j=1,3),i=1,ntri)
!
  write ( *, '(a)' ) '\\begin{figure}[htbp]'
  write ( *, '(a)' ) '\\begin{center}'
  write ( *, '(a)' ) '\\mbox{\\beginpicture'
  write ( *, '(a)' ) '\\setcoordinatesystem units <1.0cm,1.0cm>'
  write ( *,610) xmin,xmax,ymin,ymax
  write ( *, '(a)' ) '\\setlinear'

  do i = 1, ntri

   do j = 1, 3

      if ( tnbr(j,i) < i ) then

        if ( j <= 2 ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        a = til(j,i)
        b = til(jp1,i)
        if ( vcl(1,a) >= xmin .and. vcl(1,a) <= xmax .and. &
               vcl(2,a) >= ymin .and. vcl(2,a) <= ymax .and. &
               vcl(1,b) >= xmin .and. vcl(1,b) <= xmax .and. &
               vcl(2,b) >= ymin .and. vcl(2,b) <= ymax) then
              write ( *,620) vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b)
         end if

       end if

     end do

  end do

  write ( *, '(a)' ) '\\endpicture}'
  write ( *, '(a)' ) '\\end{center}'
  write ( *, '(a)' ) '\\end{figure}'

   40 continue

  610 format ('\\setplotarea x from ',f9.4,' to ',f9.4,', y from ',f9.4, &
    ' to ',f9.4)
  620 format ('\\plot ',4f15.7,' /')

  return
end subroutine



subroutine pitri
!
!***********************************************************************
!
!! PITRI generates PIC commands from DRTRI output.
!
  implicit none
!
  integer, parameter :: maxti = 8000
  integer, parameter :: maxvc = 5000
!
  integer a
  integer b
  integer i
  integer j
  integer jp1
  integer msglvl
  integer nh
  integer nhola
  integer nhole
  integer npolg
  integer ntri
  integer nvc
  integer nvcin
  integer nvert
  integer t
  integer til(3,maxti)
  integer tnbr(3,maxti)
  double precision vcl(2,maxvc)
  double precision x
  double precision xmax
  double precision xmin
  double precision ymax
  double precision ymin
!
  read ( *,*,end=40) xmin,xmax,ymin,ymax
  read ( *,*) msglvl

  if ( msglvl /= 0 ) then
    return
  end if

  read ( *, * ) nvc
  read ( *, * ) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *, * ) npolg,nhole
  read ( *, * ) (t,i=1,npolg+nhole)
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) nvert
  read ( *, * ) ((t,j=1,5),x,i=1,nvert)
  read ( *, * ) nhola,nh
  read ( *, * ) (t,i=1,nh)
  nvcin = nvc
  read ( *, * ) nvc
  read ( *, * ) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *, * ) npolg
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) nvert
  read ( *, * ) ((t,j=1,5),x,t,t,i=1,nvert)
  nvcin = nvc
  read ( *, * ) nvc
  read ( *, * ) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) ntri
  read ( *, * ) (t,(til(j,i),j=1,3),(tnbr(j,i),j=1,3),i=1,ntri)

  write ( *, '(a)' ) '.ps 5.5i'

  do i = 1, ntri

    do j = 1, 3

      if ( tnbr(j,i) < i ) then

        if ( j <= 2 ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        a = til(j,i)
        b = til(jp1,i)

        if ( vcl(1,a) >= xmin .and. vcl(1,a) <= xmax .and. &
               vcl(2,a) >= ymin .and. vcl(2,a) <= ymax .and. &
               vcl(1,b) >= xmin .and. vcl(1,b) <= xmax .and. &
               vcl(2,b) >= ymin .and. vcl(2,b) <= ymax) then
              write ( *,620) -vcl(2,a),vcl(1,a),-vcl(2,b),vcl(1,b)
        end if
!
!  Rotate by 90 degrees in above line.
!

      end if

    end do

  end do

  write ( *, '(a)' ) '.pe'

40    continue

  620 format ('line from ',f11.7,', ',f11.7,' to ', f11.7,', ',f11.7)

  return
end subroutine



subroutine test20
!
!***********************************************************************
!
!! TEST20 tests XEDGE.
!
  implicit none
!
  logical intsct
  integer mode
  double precision xu
  double precision xv1
  double precision xv2
  double precision xw1
  double precision xw2
  double precision yu
  double precision yv1
  double precision yv2
  double precision yw1
  double precision yw2
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  XEDGE determines whether two edges,'
  write ( *, '(a)' ) '  or an edge and a ray, intersect.'
  write ( *, '(a)' ) '  (An edge is a finite line segment.)'
  write ( *, '(a)' ) '  (A ray is a semi-infinite line segment.)'

  mode = 0

  xv1 = 3.0D+00
  yv1 = 0.0D+00
  xv2 = 3.0D+00
  yv2 = 2.0D+00

  xw1 = 0.0D+00
  yw1 = 0.0D+00
  xw2 = 6.0D+00
  yw2 = 2.0D+00

  write ( *, '(a)' ) ' '

  if ( mode == 0 ) then
    write ( *, '(a)' ) '  Edge 1 is from'
    write ( *, * ) '  (', xv1, ',', yv1, ') to'
    write ( *, * ) '  (', xv2, ',', yv2, ').'
  else
    write ( *, '(a)' ) '  Ray 1 is from'
    write ( *, * ) '  (', xv1, ',', yv1, ') through'
    write ( *, * ) '  (', xv2, ',', yv2, ').'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edge 2 is from'
  write ( *, * ) '  (', xw1, ',', yw1, ') to'
  write ( *, * ) '  (', xw2, ',', yw2, ').'

  call xedge ( mode, xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, xu, yu, intsct )

  if ( .not. intsct ) then
    write ( *, '(a)' ) 'The objects either do not intersect'
    write ( *, '(a)' ) 'or are parallel, or coincide.'
  else
    write ( *, '(a,2g14.6)' ) 'The objects intersect at ', xu, yu
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '(Expecting the answer (3,1) ).'

  return
end subroutine



subroutine test21
!
!***********************************************************************
!
!! TEST21 tests XLINE.
!
  implicit none
!
  double precision dv
  double precision dw
  logical parall
  double precision xu
  double precision xv1
  double precision xv2
  double precision xw1
  double precision xw2
  double precision yu
  double precision yv1
  double precision yv2
  double precision yw1
  double precision yw2
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  XLINE finds the intersection of two lines.'
  write ( *, '(a)' ) '  Each line is defined as the line a given'
  write ( *, '(a)' ) '  distance to the left of a line through two'
  write ( *, '(a)' ) '  points.'

  xv1 = 0.0D+00
  yv1 = 0.0D+00
  xv2 = 0.0D+00
  yv2 = 1.0D+00
  dv = - 6.0D+00

  xw1 = 0.0D+00
  yw1 = 0.0D+00
  xw2 = 3.0D+00
  yw2 = 1.0D+00
  dw = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, * ) '  Line 1 is ', dv, ' units left of the line'
  write ( *, * ) '  through (', xv1, ',', yv1, ') and'
  write ( *, * ) '  (', xv2, ',', yv2, ').'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, * ) '  Line 2 is ', dw, ' units left of the line'
  write ( *, * ) '  through (', xw1, ',', yw1, ') and'
  write ( *, * ) '  (', xw2, ',', yw2, ').'

  call xline ( xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, dv, dw, xu, yu, parall )

  write ( *, '(a)' ) ' '

  if ( parall ) then
    write ( *, '(a)' ) '  The lines are parallel or coincide.'
  else
    write ( *, '(a,2g14.6)' ) '  The lines intersect at ', xu, yu
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Expecting the answer (6,2) ).'

  return
end subroutine



function angle ( xa, ya, xb, yb, xc, yc )
!
!*******************************************************************************
!
!! ANGLE computes the interior angle at a vertex defined by 3 points.
!
!
!  Discussion:
!
!    ANGLE computes the interior angle, in radians, at vertex
!    (XB,YB) of the chain formed by the directed edges from
!    (XA,YA) to (XB,YB) to (XC,YC).  The interior is to the
!    left of the two directed edges.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XA, YA, XB, YB, XC, YC, the coordinates of the 
!    vertices.
!
!    Output, double precision ANGLE, the interior angle formed by
!    the vertex, in radians, between 0 and 2*PI.
!
  double precision angle
  double precision d_pi
  double precision t
  double precision x1
  double precision x2
  double precision xa
  double precision xb
  double precision xc
  double precision y1
  double precision y2
  double precision ya
  double precision yb
  double precision yc
!
  x1 = xa - xb
  y1 = ya - yb
  x2 = xc - xb
  y2 = yc - yb

  t = sqrt ( ( x1 * x1 + y1 * y1 ) * ( x2 * x2 + y2 * y2 ) )

  if ( t == 0.0D+00 ) then
    angle = d_pi ( )
    return
  end if

  t = ( x1 * x2 + y1 * y2 ) / t

  if ( t < -1.0D+00 ) then
    t = -1.0D+00
  else if ( t > 1.0D+00 ) then
    t = 1.0D+00
  end if

  angle = acos ( t )

  if ( x2 * y1 - y2 * x1 < 0.0D+00 ) then
    angle = 2.0D+00 * d_pi ( ) - angle
  end if

  return
end function

function areapg ( nvrt, xc, yc )
!
!*******************************************************************************
!
!! AREAPG computes twice the signed area of a simple polygon.
!
!
!  Modified:
!
!    13 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices on the boundary of 
!    the polygon.  N must be at least 3.
!
!    Input, double precision XC(NVRT), YC(NVRT), the X and Y coordinates 
!    of the vertices.
!
!    Output, double precision AREAPG, twice the signed area of the polygon,
!    which will be positive if the vertices were listed in counter clockwise
!    order, and negative otherwise.
!
  integer nvrt
!
  double precision areapg
  integer i
  double precision sum2
  double precision xc(nvrt)
  double precision yc(nvrt)
!
  sum2 = xc(1) * ( yc(2) - yc(nvrt) )

  do i = 2, nvrt-1
    sum2 = sum2 + xc(i) * ( yc(i+1) - yc(i-1) )
  end do

  sum2 = sum2 + xc(nvrt) * ( yc(1) - yc(nvrt-1) )

  areapg = sum2

  return
end function

function areatr ( xa, ya, xb, yb, xc, yc )
!
!*******************************************************************************
!
!! AREATR computes twice the signed area of a triangle.
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XA, YA, XB, YB, XC, YC, the coordinates of the 
!    vertices.
!
!    Output, double precision AREATR, twice the signed area of the triangle.
!    This will be positive if the vertices are listed in counter clockwise 
!    order.
!
  double precision areatr 
  double precision xa
  double precision xb
  double precision xc
  double precision ya
  double precision yb
  double precision yc
!
  areatr = ( xb - xa ) * ( yc - ya ) - ( xc - xa ) * ( yb - ya )

  return
end function


subroutine bedgmv ( nvc, npolg, nvert, maxvc, h, vcl, hvl, pvl, vstart, vnum, &
  ierror )
!
!*******************************************************************************
!
!! BEDGMV generates boundary edge mesh vertices.
!
!
!  Purpose: 
!
!    Generate mesh vertices on boundary of convex polygons
!    of decomposition with spacing determined by H array.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input/output, integer NVC, the number of coordinates or positions used 
!    in VCL array.
!
!    Input, integer NPOLG, the number of polygons or positions used in HVL array.
!
!    Input, integer NVERT, the number of vertices or positions used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, double precision H(1:NPOLG), the spacing of mesh vertices for 
!    convex polygons.
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input, integer HVL(1:NPOLG, the head vertex list.
!
!    Input, integer PVL(1:4,1:NVERT), the polygon vertex list.
!
!    Output, integer VSTART(1:NVERT), the start location in VCL for mesh 
!    vertices on each edge in PVL if there are any, else 0.
!
!    Output, integer VNUM(1:NVERT), the number of mesh vertices on interior 
!    of each edge in PVL; entry is negated if mesh vertices are listed in 
!    backward order in VCL.
!
!    Output, integer IERROR, is set to 3 on error.
!
  integer maxvc
  integer npolg
  integer nvert
!
  double precision dx
  double precision dy
  integer, parameter :: edgv = 4
  double precision h(npolg)
  double precision hh
  integer hvl(npolg)
  integer i
  integer ia
  integer ierror
  integer j
  integer k
  integer l
  double precision leng
  integer, parameter :: loc = 1
  integer m
  integer nvc
  integer, parameter :: polg = 2
  integer pvl(4,nvert)
  integer, parameter :: succ = 3
  integer u
  integer v
  integer vstart(nvert)
  integer vnum(nvert)
  double precision vcl(2,maxvc)
  double precision x
  double precision y
!
  ierror = 0
  vstart(1:nvert) = -1

  do k = 1, npolg

    i = hvl(k)

    do

      j = pvl(succ,i)

      if ( vstart(i) == -1 ) then

        u = pvl(loc,i)
        v = pvl(loc,j)
        x = vcl(1,u)
        y = vcl(2,u)
        leng = sqrt ( ( vcl(1,v) - x )**2 + ( vcl(2,v) - y )**2 )
        ia = pvl(edgv,i)

        if ( ia <= 0 ) then
          hh = h(k)
        else
          hh = sqrt ( h(k) * h(pvl(polg,ia)) )
        end if

        if ( hh == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BEDGMV - Fatal error!'
          write ( *, '(a)' ) '  HH = 0.'
          stop
        end if

        l = int ( leng / hh )

        if ( ( leng / hh ) - dble ( l ) > dble(l) / dble ( 2 * l + 1 ) ) then
          l = l + 1
        end if

        if ( l <= 1 ) then

          vstart(i) = 0
          vnum(i) = 0

        else

          dx = ( vcl(1,v) - x ) / dble ( l )
          dy = ( vcl(2,v) - y ) / dble ( l )
          l = l - 1

          if ( nvc + l > maxvc ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'BEDGMV - Fatal error!'
            write ( *, '(a)' ) '  NVC+L > MAXVC.'
            ierror = 3
            return
          end if

          vstart(i) = nvc + 1
          vnum(i) = l

          do m = 1, l
            x = x + dx
            y = y + dy
            nvc = nvc + 1
            vcl(1,nvc) = x
            vcl(2,nvc) = y
          end do

        end if

        if ( ia > 0 ) then
          vstart(ia) = vstart(i)
          vnum(ia) = -vnum(i)
        end if

      end if

      i = j

      if ( i == hvl(k) ) then
        exit
      end if

    end do

  end do

  return
end subroutine



subroutine bnsrt2 ( binexp, n, a, map, bin, iwk )
!
!*******************************************************************************
!
!! BNSRT2 bin sorts N points in 2D into increasing bin order.
!
!
!  Purpose: 
!
!    Use a bin sort to obtain the permutation of N 2-dimensional
!    double precision points so that points are in increasing bin
!    order, where the N points are assigned to about N**BINEXP bins.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer BINEXP, the exponent for the number of bins.
!
!    Input, integer N, the number of points.
!
!    Input, double precision A(2,*), the points to be binned.
!
!    Input/output, integer MAP(N); on input, the points of A with indices 
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, MAP has
!    been permuted so bin of MAP(1) <= bin of MAP(2) <= ... <= bin of MAP(N).
!
!    Workspace, integer BIN(N), used for bin numbers and permutation of 1 to N.
!
!    Workspace, intger IWK(N), used for copy of MAP array.
!
  integer n
!
  double precision a(2,*)
  integer bin(n)
  double precision binexp
  double precision dx
  double precision dy
  integer i
  integer iwk(n)
  integer j
  integer k
  integer l
  integer map(n)
  integer nside
  double precision xmax
  double precision xmin
  double precision ymax
  double precision ymin
!
  nside = int ( dble ( n )**( binexp / 2.0D+00 ) + 0.5D+00 )

  if ( nside <= 1 ) then
    return
  end if

  xmin = a(1,map(1))
  ymin = a(2,map(1))
  xmax = xmin
  ymax = ymin
  do i = 2, n
    j = map(i)
    xmin = min ( xmin, a(1,j) )
    xmax = max ( xmax, a(1,j) )
    ymin = min ( ymin, a(2,j) )
    ymax = max ( ymax, a(2,j) )
  end do

  dx = 1.0001D+00 * ( xmax - xmin ) / dble ( nside )
  dy = 1.0001D+00 * ( ymax - ymin ) / dble ( nside )

  if ( dx == 0.0D+00 ) then
    dx = 1.0D+00
  end if

  if ( dy == 0.0D+00 ) then
    dy = 1.0D+00
  end if

  do i = 1, n
    j = map(i)
    iwk(i) = j
    map(i) = i
    k = int ( ( a(1,j) - xmin ) / dx )
    l = int ( ( a(2,j) - ymin ) / dy )
    if ( mod ( k, 2 ) == 0 ) then
      bin(i) = k * nside + l
    else
      bin(i) = ( k + 1 ) * nside - l - 1
    end if
  end do

  call ihpsrt ( 1, n, 1, bin, map )

  bin(1:n) = map(1:n)

  do i = 1, n
    map(i) = iwk(bin(i))
  end do

  return
end subroutine



function cmcirc ( x0, y0, x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! CMCIRC determines whether a point lies within a circle through 3 points.
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X0, Y0, the coordinates of the point to
!    be tested.
!
!    Input, double precision X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    three points that define a circle.
!
!    Output, integer CMCIRC, reports the test results:
!     2, if the three vertices are collinear,
!     1, if (X0,Y0) is inside the circle,
!     0, if (X0,Y0) is on the circle,
!    -1, if (X0,Y0) is outside the circle.
!
  double precision a11
  double precision a12
  double precision a21
  double precision a22
  double precision b1
  double precision b2
  integer cmcirc
  double precision det
  double precision diff
  double precision rsq
  double precision tol
  double precision tolabs
  double precision xc
  double precision yc
  double precision x0
  double precision x1
  double precision x2
  double precision x3
  double precision y0
  double precision y1
  double precision y2
  double precision y3
!
  tol = 100.0D+00 * epsilon ( tol )
  cmcirc = 2
  a11 = x2 - x1
  a12 = y2 - y1
  a21 = x3 - x1
  a22 = y3 - y1
  tolabs = tol * max ( abs ( a11), abs ( a12), abs ( a21), abs ( a22) )
  det = a11 * a22 - a21 * a12

  if ( abs ( det ) <= tolabs ) then
    return
  end if

  b1 = a11 * a11 + a12 * a12
  b2 = a21 * a21 + a22 * a22
  det = 2.0D+00 * det
  xc = ( b1 * a22 - b2 * a12 ) / det
  yc = ( b2 * a11 - b1 * a21 ) / det
  rsq = xc * xc + yc * yc
  diff = ( ( x0 - x1 - xc)**2 + ( y0 - y1 - yc )**2 ) - rsq
  tolabs = tol * rsq

  if ( diff < - tolabs ) then
    cmcirc = 1
  else if ( diff > tolabs ) then
    cmcirc = -1
  else
    cmcirc = 0
  end if

  return
end function



subroutine cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, &
  maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )
!
!*******************************************************************************
!
!! CVDEC2 decomposes a polygonal region into convex polygons.
!
!
!  Purpose: 
!
!    Decompose general polygonal region (which is decomposed
!    into simple polygons on input) into convex polygons using
!    vertex coordinate list, head vertex list, and polygon vertex
!    list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGSPC, the angle spacing parameter in radians 
!    used in controlling vertices to be considered as an endpoint of a separator.
!
!    Input, double precision ANGTOL, the angle tolerance parameter in radians 
!    used in accepting separator(s).
!
!    Input/output, integer NVC, the number of vertex coordinates or positions 
!    used in VCL.
!
!    Input/output, integer NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array, should be >=
!    number of vertex coordinates required for decomposition.
!
!    Input, integer MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be >= number of polygons required for decomposition.
!
!    Input, integer MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be >= number of polygon vertices required for decomposition.
!
!    Input, integer MAXIW, the maximum size available for IWK array; should be 
!    about 3 times maximum number of vertices in any polygon.
!
!    Input, integer MAXWK, the maximum size available for WK array; should be 
!    about 5 times maximum number of vertices in any polygon.
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer REGNUM(1:NPOLG), region numbers.
!
!    Input/output, integer HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT),
!    the polygon vertex list and interior angles; see routine DSPGDC for 
!    more details.  Note that the data structures should be as output from 
!    routine SPDEC2.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  For abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 206, 207, 208, 209, 210, or 212.
!
  integer maxhv
  integer maxiw
  integer maxpv
  integer maxvc
  integer maxwk
!
  double precision angspc
  double precision angtol
  integer hvl(maxhv)
  double precision iang(maxpv)
  integer ierror
  integer iwk(maxiw)
  integer npolg
  integer nvc
  integer nvert
  double precision d_pi
  double precision piptol
  integer pvl(4,maxpv)
  integer regnum(maxhv)
  double precision tol
  integer v
  double precision vcl(2,maxvc)
  integer w1
  integer w2
  double precision wk(maxwk)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  For each reflex vertex, resolve it with one or two separators
!  and update VCL, HVL, PVL, IANG.
!
  piptol = d_pi ( ) + tol
  v = 1

  do

    if ( v > nvert ) then
      exit
    end if

    if ( iang(v) > piptol ) then

      call resvrt ( v, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
        maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call insed2 ( v ,w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
        pvl, iang, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( w2 > 0 ) then
        call insed2 ( v, w2, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
          pvl, iang, ierror )
      end if

      if ( ierror /= 0 ) then
        return
      end if

    end if

    v = v + 1

  end do

  return
end subroutine



subroutine cvdtri ( inter, ldv, nt, vcl, til, tedg, sptr, ierror )
!
!*******************************************************************************
!
!! CVDTRI converts boundary triangles to Delaunay triangles.
!
!
!  Purpose: 
!
!    Convert triangles in strip near boundary of polygon
!    or inside polygon to Delaunay triangles.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical INTER, is .TRUE. if and only if there is at least 
!    one interior mesh vertex.
!
!    Input, integer LDV, the leading dimension of VCL in calling routine.
!
!    Input, integer NT, the number of triangles in strip or polygon.
!
!    Input, VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer TIL(1:3,1:NT), the triangle incidence list.
!
!    Input/output, integer TEDG(1:3,1:NT) - TEDG(J,I) refers to edge with 
!    vertices TIL(J:J+1,I) and contains index of merge edge or
!    > NT for edge of chains.
!
!    Workspace, SPTR(1:NT) - SPTR(I) = -1 if merge edge I is not in LOP stack,
!    else >= 0 and pointer (index of SPTR) to next edge in
!    stack (0 indicates bottom of stack).
!
!    Output, integer IERROR, error flag.  On abnormal return:
!    IERROR is set to 231.
!
  integer ldv
  integer nt
!
  integer e
  integer ierror
  integer ind(2)
  logical inter
  integer itr(2)
  integer k
  integer mxtr
  logical sflag
  integer sptr(nt)
  integer tedg(3,nt)
  integer til(3,nt)
  integer top
  double precision vcl(ldv,*)
!
  ierror = 0
  sflag = .true.
  sptr(1:nt) = -1

  do k = 1, nt

    mxtr = k + 1

    if ( k == nt ) then
      if ( .not. inter ) then
        return
      end if
      mxtr = nt
      sflag = .false.
    end if

    top = k
    sptr(k) = 0

    do

      e = top
      top = sptr(e)

      call fndtri ( e, mxtr, sflag, tedg, itr, ind, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call lop ( itr, ind, k, top, ldv, vcl, til, tedg, sptr )

      if ( top <= 0 ) then
        exit
      end if

    end do

  end do

  return
end subroutine



function degrees_to_radians ( angle )
!
!*******************************************************************************
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision ANGLE, an angle in degrees.
!
!    Output, double precision DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  double precision angle
  double precision degrees_to_radians
  double precision d_pi
!
  degrees_to_radians = ( angle / 180.0D+00 ) * d_pi ( )

  return
end function



subroutine delaunay_print ( num_pts, xc, num_tri, nodtri, tnbr )
!
!*******************************************************************************
!
!! DELAUNAY_PRINT prints out information defining a Delaunay triangulation.
!
!
!  Modified:
!
!    08 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NUM_PTS, the number of points.
!
!    Input, double precision XC(2,NUM_PTS), the point coordinates.
!
!    Input, integer NUM_TRI, the number of triangles.
!
!    Input, integer NODTRI(3,NUM_TRI), the nodes that make up the triangles.
!
!    Input, integer TNBR(3,NUM_TRI), the triangle neighbors on each side.
!
  integer num_pts
  integer num_tri
!
  integer i
  integer i_wrap
  integer j
  integer k
  integer n1
  integer n2
  integer nodtri(3,num_tri)
  integer s
  integer t
  integer tnbr(3,num_tri)
  double precision xc(2,num_pts)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DELAUNAY_PRINT'
  write ( *, '(a)' ) '  Information defining a Delaunay triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points is ', num_pts

  call dmat_print ( num_pts, num_pts, 2, transpose ( xc ), &
    '  Point coordinates (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of triangles is ', num_tri
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three points are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the points'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call imat_print ( num_tri, num_tri, 3, transpose ( nodtri ), &
    '  Nodes that make up triangles (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call imat_print ( num_tri, num_tri, 3, transpose ( tnbr ), &
    '  Indices of neighboring triangles (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of boundary points (and segments) is ', &
    2 * num_pts - num_tri - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  # Tri Side  N1  N2'
  write ( *, '(a)' ) ' '
  k = 0

  do i = 1, num_tri
    do j = 1, 3
      if ( tnbr(j,i) < 0 ) then
        s = - tnbr(j,i)
        t = s / 3
        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = nodtri(s,t)
        n2 = nodtri(i_wrap(s+1,1,3),t)
        write ( *, '(5i4)' ) k, t, s, n1, n2
      end if
    end do
  end do

  return
end subroutine


subroutine dhpsrt ( k, n, lda, a, map )
!
!*******************************************************************************
!
!! DHPSRT sorts points into lexicographic order using heap sort
!
!
!  Discussion:
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    points so that the points are in lexicographic increasing order.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Modified:
!
!    19 February 2001
!
!  Parameters:
!
!    Input, integer K, the dimension of the points (for instance, 2
!    for points in the plane).
!
!    Input, integer N, the number of points.
!
!    Input, integer LDA, the leading dimension of array A in the calling
!    routine; LDA should be at least K.
!
!    Input, double precision A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted.
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  integer lda
  integer n
!
  double precision a(lda,*)
  integer i
  integer k
  integer map(n)
!
  do i = n/2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    call i_swap ( map(1), map(i) )
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end subroutine


function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! DIAEDG chooses one of the diagonals of a quadrilateral.
!
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge 
!    that should be chosen, based on the circumcircle criterion, where 
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple 
!    quadrilateral in counterclockwise order.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  double precision ca
  double precision cb
  integer diaedg
  double precision dx10
  double precision dx12
  double precision dx30
  double precision dx32
  double precision dy10
  double precision dy12
  double precision dy30
  double precision dy32
  double precision s
  double precision tol
  double precision tola
  double precision tolb
  double precision x0
  double precision x1
  double precision x2
  double precision x3
  double precision y0
  double precision y1
  double precision y2
  double precision y3
!
  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( ca > tola .and. cb  >  tolb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( s > tola ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end function


subroutine diam2 ( nvrt, xc, yc, i1, i2, diamsq, ierror )
!
!*******************************************************************************
!
!! DIAM2 finds the diameter of a convex polygon.
!
!
!  Purpose: 
!
!    Find the diameter of a convex polygon with vertices
!    given in counter clockwise order and with all interior angles < PI.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices.
!
!    Input, double precision XC(NVRT),YC(NVRT), the vertex coordinates in 
!    counter-clockwise order.
!
!    Output, integer I1, I2 , indices in XC, YC of the diameter edge; the
!    diameter is from (XC(I1),YC(I1)) to (XC(I2),YC(I2)).
!
!    Output, double precision DIAMSQ, the square of the diameter.
!
!    Output, integer IERROR, an error flag.
!    0, no error was detected.
!    200, an error was detected.
!
  integer nvrt
!
  double precision area1
  double precision area2
  double precision areatr
  double precision c1mtol
  double precision c1ptol
  double precision diamsq
  double precision dist
  integer i1
  integer i2
  integer ierror
  integer j
  integer jp1
  integer k
  integer kp1
  integer m
  double precision tol
  double precision xc(nvrt)
  double precision yc(nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Find the first vertex which is farthest from the edge connecting
!  vertices with indices NVRT, 1.
!
  c1mtol = 1.0D+00 - tol
  c1ptol = 1.0D+00 + tol
  j = nvrt
  jp1 = 1
  k = 2
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

    if ( area2 <= area1 * c1ptol ) then
      exit
    end if

    area1 = area2
    k = k + 1

  end do

  m = k
  diamsq = 0.0D+00
!
!  Find diameter = maximum distance of antipodal pairs.
!
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    kp1 = k + 1
    if ( kp1 > nvrt) kp1 = 1

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(kp1), yc(kp1) )

    if ( area2 > area1 * c1ptol ) then
      k = k + 1
      area1 = area2
    else if ( area2 < area1 * c1mtol ) then
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
    else
      k = k + 1
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
    end if

    if ( j > m .or. k > nvrt ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIAM2 - Fatal error!'
      write ( *, '(a)' ) '  J < M or K > NVRT.'
      ierror = 200
      return
    end if

    dist = ( xc(j) - xc(k) )**2 + ( yc(j) - yc(k) )**2

    if ( dist > diamsq ) then
      diamsq = dist
      i1 = j
      i2 = k
    end if

    if ( j == m .and. k == nvrt ) then
      exit
    end if

  end do

  return
end subroutine



function dless ( k, p, q )
!
!*******************************************************************************
!
!! DLESS determine whether P is lexicographically less than Q.
!
!
!  Discussion:
!
!    P and Q are K-dimensional points.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer K, the spatial dimension of the points.
!
!    Input, double precision P(K), Q(K), the points to be compared.
!
!    Output, logical RLESS, is TRUE if P < Q, FALSE otherwise.
!
  integer k
!
  double precision cmax
  logical dless
  integer i
  double precision p(k)
  double precision q(k)
  double precision tol
!
  tol = 100.0D+00 * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( abs ( p(i) - q(i) ) > tol * cmax .and. cmax > tol ) then

      if ( p(i) < q(i) ) then
        dless = .true.
      else
        dless = .false.
      end if

      return
    end if

  end do

  dless = .false.

  return
end function


subroutine dmat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! DMAT_PRINT prints a double precision matrix.
!
!
!  Modified:
!
!    08 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, double precision A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer lda
  integer n
!
  double precision a(lda,n)
  integer i
  integer j
  integer jhi
  integer jlo
  integer m
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end subroutine



subroutine dsftdw ( l, u, k, lda, a, map )
!
!*******************************************************************************
!
!! DSFTDW sifts A(*,MAP(L)) down a heap of size U.
!
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer L, U, the lower and upper indices of part of the heap.
!
!    Input, integer K, the spatial dimension of the points.
!
!    Input, integer LDA, the leading dimension of A in the calling routine.
!
!    Input, double precision A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted.
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  integer lda
!
  double precision a(lda,*)
  logical dless
  integer i
  integer j
  integer k
  integer l
  integer map(*)
  integer t
  integer u
!
  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end subroutine



subroutine dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxho, nvc, npolg, &
  nvert, nhola, regnum, hvl, pvl, iang, holv, ierror )
!
!*******************************************************************************
!
!! DSMCPR initializes the polygonal decomposition data structure.
!
!
!  Purpose: 
!
!    Initialize the polygonal decomposition data structure
!    given a multiply-connected polygonal region with 1 outer
!    boundary curve and 0 or more inner boundary curves of holes.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NHOLE, the number of holes in region.
!
!    Input, integer NVBC(1:NHOLE+1), the number of vertices per boundary curve; 
!    first boundary curve is the outer boundary of the region.
!
!    Input, double precision VCL(1:2,1:NVC), the vertex coordinates of boundary 
!    curves in counter clockwise order; NVC = NVBC(1) + ... + NVBC(NHOLE+1); 
!    positions 1 to NVBC(1) of VCL contain the vertex coordinates of the
!    outer boundary in counter clockwise order; positions NVBC(1)+1 to
!    NVBC(1)+NVBC(2) contain the vertex coordinates of the
!    first hole boundary in counter clockwise order, etc.
!
!    Input, integer MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be >= NHOLE + 1.
!
!    Input, integer MAXPV, the maximum size available for PVL, IANG arrays;
!    should be >= NVC.
!
!    Input, integer MAXHO, the maximum size available for HOLV array; should be
!    >= NHOLE*2.
!
!    Output, integer NVC, the number of vertex coordinates, set to sum of 
!    NVBC(I).
!
!    Output, integer NPOLG, the number of polygonal subregions, set to 1.
!    For consistency with DSPGDC.
!
!    Output, integer NVERT, the number of vertices in PVL, set to NVC.
!    For consistency with DSPGDC.
!
!    Output, integer NHOLA, number of attached holes, set to 0.
!    For consistency with DSPGDC.
!
!    Output, integer REGNUM(1:1), region number of only subregion, set to 1
!    For consistency with DSPGDC.
!
!    Output, integer HVL(1:NHOLE+1), the head vertex list; first entry is the 
!    head vertex (index in PVL) of outer boundary curve; next
!    NHOLE entries contain the head vertex of a hole.
!
!    Output, integer PVL(1:4,1:NVC), IANG(1:NVC), the polygon vertex list and 
!    interior angles; vertices of outer boundary curve are in counter clockwise 
!    order followed by vertices of each hole in CW hole; vertices
!    of each polygon are in a circular linked list; see
!    routine DSPGDC for more details of this data structure.
!
!    Output, integer HOLV(1:NHOLE*2), the indices in PVL of top and bottom 
!    vertices of holes; first (last) NHOLE entries are for top (bottom)
!    vertices; top (bottom) vertices are sorted in decreasing
!    (increasing) lexicographic (y,x) order of coordinates.
!
!    Output, integer IERROR, error flag.
!    For abnormal return, IERROR is set to 2, 4, or 5.
!
  integer maxho
  integer maxpv
  integer nhole
!
  double precision angle
  integer, parameter :: edgv = 4
  integer i
  double precision iang(maxpv)
  integer ierror
  integer iv
  integer ivs
  integer j
  integer, parameter :: loc = 1
  integer lv
  integer lvp
  integer lvs
  integer maxhv
  integer nhola
  integer npolg
  integer nv
  integer nvc
  integer nvert
  integer nvs
  integer hvl(nhole+1)
  integer holv(maxho)
  integer nvbc(nhole+1)
  integer, parameter :: polg = 2
  integer pvl(4,maxpv)
  integer regnum(1)
  integer, parameter :: succ = 3
  double precision vcl(2,*)
!
  ierror = 0
  nvc = sum ( nvbc(1:nhole+1) )
  npolg = 1
  nvert = nvc
  nhola = 0
  regnum(1) = 1

  if ( nhole + 1 > maxhv ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSMCPR - Fatal error!'
    write ( *, '(a)' ) '  NHOLE+1 > MAXHV.'
    ierror = 4
    return
  end if

  if ( nvc > maxpv ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSMCPR - Fatal error!'
    write ( *, '(a)' ) '  NVC > MAXPV.'
    ierror = 5
    return
  end if

  if ( nhole + nhole > maxho ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSMCPR - Fatal error!'
    write ( *, '(a)' ) '  2 * NHOLE > MAXHO.'
    ierror = 2
    return
  end if
!
!  Initialize the HVL and PVL arrays.
!
  hvl(1) = 1
  nv = nvbc(1)

  do i = 1, nv
    pvl(loc,i) = i
    pvl(polg,i) = 1
    pvl(succ,i) = i + 1
    pvl(edgv,i) = 0
  end do

  pvl(succ,nv) = 1

  do j = 1, nhole
    hvl(j+1) = nv + 1
    nvs = nv + nvbc(j+1)
    do i = nv+1, nvs
      pvl(loc,i) = i
      pvl(polg,i) = 1
      pvl(succ,i) = i - 1
      pvl(edgv,i) = 0
    end do
    pvl(succ,nv+1) = nvs
    nv = nvs
  end do
!
!  Initialize the IANG array.
!
  do i = 1, nhole+1

    j = hvl(i)
    lvp = pvl(loc,j)
    iv = pvl(succ,j)
    lv = pvl(loc,iv)

    do

      ivs = pvl(succ,iv)
      lvs = pvl(loc,ivs)
      iang(iv) = angle ( vcl(1,lvp), vcl(2,lvp), vcl(1,lv), vcl(2,lv), &
        vcl(1,lvs), vcl(2,lvs) )

      if ( iv == j ) then
        exit
      end if

      lvp = lv
      iv = ivs
      lv = lvs

    end do

  end do
!
!  Initialize HOLV array.
!
  if ( nhole > 0 ) then
    call holvrt ( nhole, vcl, hvl(2), pvl, holv )
  end if

  return
end subroutine


subroutine dsmdf2 ( hflag, nvc, npolg, maxwk, vcl, hvl, pvl, iang, ivrt, &
  xivrt, widsq, edgval, vrtval, area, wk, ierror )
!
!*******************************************************************************
!
!! DSMDF2 sets up a data structure for a heuristic mesh distribution.
!
!
!  Purpose: 
!
!    Set up the data structure for heuristic mesh distribution
!    function from data structure for convex polygon decomposition
!    if HFLAG is .TRUE., else set up only IVRT and XIVRT.
!
!    Also compute areas of convex polygons.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical HFLAG, set to .TRUE. if data structure is to be constructed,
!    .FALSE. if only IVRT, XIVRT, AREA are to be computed.
!
!    Input, integer NVC, the number of vertex coordinates in VCL array.
!
!    Input, integer NPOLG, the number of polygonal subregions in HVL array.
!
!    Input, integer MAXWK, the maximum size available for WK array; should be
!    2 times maximum number of vertices in any polygon.
!
!    Input, VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input, integer HVL(1:NPOLG), the head vertex list.
!
!    Input, integer PVL(1:4,1:*), double precision IANG(1:*), the polygon vertex
!    list, interior angles.
!
!    Output, integer IVRT(1:*), the indices of polygon vertices in VCL, ordered 
!    by polygon; same size as PVL.  For heuristic MDF data structure.
!
!    Output, XIVRT(1:NPOLG+1), the pointer to first vertex of each polygon
!    in IVRT; vertices of polygon K are IVRT(I) for I from
!    XIVRT(K) to XIVRT(K+1)-1.  For heuristic MDF data structure.
!
!    Output, double precision WIDSQ(1:NPOLG), the square of width of convex 
!    polygons.  For heuristic MDF data structure.
!
!    Output, double precision EDGVAL(1:*), the value associated with each 
!    edge of decomposition; same size as PVL.  For heuristic MDF data structure.
!
!    Output, double precision VRTVAL(1:NVC), the value associated with each 
!    vertex of decomposition.  For heuristic MDF data structure.
!
!    Output, double precision AREA(1:NPOLG), the area of convex polygons.
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 7 or 201.
!
  integer maxwk
  integer npolg
  integer nvc
!
  double precision area(npolg)
  double precision areapg
  integer, parameter :: edgv = 4
  double precision edgval(*)
  logical hflag
  integer hvl(npolg)
  integer i
  double precision iang(*)
  integer ierror
  integer il
  integer ivrt(*)
  integer j
  integer jl
  integer k
  integer l
  integer, parameter :: loc = 1
  integer m
  integer nvrt
  double precision d_pi
  double precision pimtol
  integer, parameter :: polg = 2
  integer pvl(4,*)
  double precision s
  integer, parameter :: succ = 3
  double precision tol
  double precision vcl(2,nvc)
  double precision vrtval(nvc)
  double precision widsq(npolg)
  double precision wk(maxwk)
  integer xc
  integer xivrt(npolg+1)
  integer yc
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Compute area and square of width of polygons.
!
  pimtol = d_pi ( ) - tol

  do k = 1, npolg

    nvrt = 0
    i = hvl(k)

    do

      if ( iang(i) < pimtol ) then
        nvrt = nvrt + 1
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    if ( nvrt + nvrt > maxwk ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DSMDF2 - Fatal error!'
      write ( *, '(a)' ) '  NVRT + NVRT > MAXWK.'
      write ( *, '(a,i6)' ) '  NVRT = ', nvrt
      write ( *, '(a,i6)' ) '  MAXWK = ', maxwk
      ierror = 7
      return
    end if

    xc = 0

    do

      if ( iang(i) < pimtol ) then
        j = pvl(loc,i)
        xc = xc + 1
        wk(xc) = vcl(1,j)
        wk(xc+nvrt) = vcl(2,j)
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    xc = 1
    yc = xc + nvrt
    area(k) = areapg ( nvrt, wk(xc), wk(yc) ) * 0.5D+00

    if ( hflag ) then

      call width2 ( nvrt, wk(xc), wk(yc), i, j, widsq(k), ierror )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DSMDF2 - Fatal error!'
        write ( *, '(a)' ) '  WIDTH2 returns an error condition.'
        return
      end if

    end if

  end do
!
!  Set up IVRT, XIVRT, EDGVAL, VRTVAL arrays.
!
  l = 1

  do k = 1, npolg

    xivrt(k) = l
    i = hvl(k)
    il = pvl(loc,i)

    do

      ivrt(l) = il
      j = pvl(succ,i)
      jl = pvl(loc,j)

      if ( hflag ) then

        s = min ( (vcl(1,jl) - vcl(1,il) )**2 + ( vcl(2,jl) - vcl(2,il) )**2, &
          widsq(k) )

        m = pvl(edgv,i)
        if ( m > 0 ) then
          s = min ( s, widsq(pvl(polg,m) ) )
        end if

        edgval(l) = s

      end if

      l = l + 1
      i = j
      il = jl

      if ( i == hvl(k) ) then
        exit
      end if

    end do

  end do

  xivrt(npolg+1) = l

  if ( .not. hflag ) then
    return
  end if

  vrtval(1:nvc) = 0.0D+00

  do k = 1, npolg

    j = xivrt(k+1) - 1
    l = j

    do i = xivrt(k),l

      il = ivrt(i)

      if ( vrtval(il) == 0.0D+00 ) then
        vrtval(il) = min ( edgval(i), edgval(j) )
      else
        vrtval(il) = min ( vrtval(il), edgval(i), edgval(j) )
      end if

      j = i

    end do

  end do

  return
end subroutine



subroutine dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv,  &
  maxho, npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, holv, htsiz,  &
  maxedg, ht, edge, map, ierror )
!
!*******************************************************************************
!
!! DSPGDC initializes the polygonal decomposition data structure.
!
!
!  Purpose: 
!
!    Initialize the polygonal decomposition data structure
!    given an initial decomposition of a polygonal region which
!    may have holes and/or cut, separator, and hole interfaces.
!    Holes and hole interfaces must be simple polygons.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVC, the number of distinct vertex coordinates in region.
!
!    Input, double precision VCL(1:2,1:NVC), the vertex coordinates of boundary 
!    curves in arbitrary order.
!
!    Input, integer INCR,  a positive integer >= NVC, e.g. 10000, added to some
!    elements of IVRT array.
!
!    Input, integer NCUR, the number of boundary curves (includes outer boundary
!    curves of subregions and boundary curves of holes
!    and hole interfaces).
!
!    Input, integer NVBC(1:NCUR), the number of vertices per boundary curve.
!
!    Input, integer ICUR(1:NCUR), indicates type and location of the curves:
!    ICUR(I) = 0 if Ith curve is outer boundary curve,
!    ICUR(I) = K if Ith curve is a hole and is inside
!      the subregion to the left of Kth curve,
!    ICUR(I) = -K if Ith curve is a hole interface and is
!      inside the subregion to the left of Kth curve.
!    K must be the index of an outer or hole interface
!    boundary curve (hole interfaces may be nested).
!    If the Ith curve is inside more than one subregion
!    due to nesting of hole interfaces, then the subregion
!    to the left of Kth curve must be the smallest
!    subregion containing the Ith curve.
!
!    Input, integer IVRT(1:NV), indices in VCL of vertices of boundary curves;
!    NV = NVBC(1) + ... + NVBC(NCUR); the vertices of each
!    boundary curve must be in counter clockwise order; the first NVBC(1)
!    positions of IVRT are used for the first curve; the
!    next NVBC(2) positions are used for second curve, etc.
!    If the Ith curve is the outer boundary of a subregion
!    determined from cut and separator interfaces, then the
!    elements of IVRT which correspond to this curve are used
!    both for an index in VCL and indicating the type of the
!    edge joining a vertex and its successor as follows.
!    Let J be in range of positions used for the Ith curve
!    and K be the index in VCL of the coordinates of a vertex
!    of the Ith curve. Consider the edge originating from this
!    vertex. IVRT(J) = -K if the edge is part of a cut or
!    separator interface (i.e. there is a subregion to right
!    of edge). IVRT(J) = K if the edge is part of the outer
!    boundary of the region (i.e. the unbounded exterior of
!    the region is to the right of edge). IVRT(J) = K + INCR
!    if the edge is part of the boundary of a hole (i.e.
!    there is a bounded area to the right of edge which is
!    not in the region. If the Ith curve is the boundary of
!    a hole or hole interface, then only IVRT(J) = K is used.
!
!    Input, integer MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be >= NCUR + (number of hole interfaces).
!
!    Input, integer MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be >= NVERT (see below).
!
!    Input, integer MAXHO, the maximum size available for HOLV array; should be
!    >= NHOLE*2 + NHOLA (see below).
!
!    Input, integer HTSIZ, the size of hash table HT; should be a prime number 
!    which is about NSC/2 where NSC is number of separator and cut
!    interface edges.
!
!    Input, integer MAXEDG, the maximum size available for EDGE array; should 
!    be at least NSC.
!
!    Output, integer NPOLG, the number of polygonal subregions, set to number 
!    of outer subregion boundaries plus number of hole interfaces.
!
!    Output, integer NVERT, the number of vertices in PVL, set to NV plus number
!    of vertices in holes and hole interfaces (< 2*NV).
!
!    Output, integer NHOLE, the number of holes and hole interfaces.
!
!    Output, integer NHOLA, the number of 'attached' holes; these holes are 
!    attached to the outer boundary of a subregion through vertices
!    or cut interfaces and have their edges in consecutive
!    order on the boundary (<= NV/4).
!
!    Output, integer REGNUM(1:NPOLG). region numbers to left of outer and hole
!    interface boundary curves, which are set to the indices
!    of ICUR or NVBC; this array may be useful in some
!    applications for identifying which original region a
!    subpolygon belongs to.
!
!    Output, HVL(1:NPOLG+NHOLE), the head vertex list; the first NPOLG
!    positions contain the head vertex (index in PVL) of an
!    outer or hole interface boundary curve in which the
!    vertices of the curve are in counter clockwise order in PVL; next
!    NHOLE positions contain the head vertex of a hole or
!    hole interface in which vertices are in CW order in PVL.
!
!    Output, PVL(1:4,1:NVERT), double precision IANG(1:NVERT), the polygon 
!    vertex list and interior angles; contains the 5 'arrays' LOC, POLG, SUCC
!    EDGV, IANG (the first 4 are integer arrays, the last
!    is a double precision array); the vertices of each
!    polygon (except for holes) are stored in counter clockwise order in a
!    circular linked list. PVL(LOC,V) is the location in VCL
!    of the coordinates of 'vertex' (index) V. IANG(V) is
!    the interior angle at vertex V. PVL(POLG,V) is polygon
!    number (index of HVL) of subregion containing vertex V
!    (this entry is different from the polygon index only
!    for holes). PVL(SUCC,V) is index in PVL of successor
!    vertex of vertex V. PVL(EDGV,V) gives information about
!    the edge joining vertices V and its successor - if the
!    edge is part of 1 polygon then PVL(EDGV,V) = 0; if the
!    edge is common to 2 polygons then PVL(EDGV,V) > 0 and
!    is equal to the index in PVL of the successor vertex
!    as represented in the other polygon; i.e. in latter
!    case, PVL(LOC,PVL(EDGV,V)) = PVL(LOC,PVL(SUCC,V)) and
!    PVL(EDGV,PVL(EDGV,V)) = V.
!
!    Output, integer HOLV(1:NHOLE*2+NHOLA), indices in PVL of top or bottom 
!    vertex of holes; first (next) NHOLE entries are for top (bottom)
!    vertices of holes and hole interfaces, with top (bottom)
!    vertices sorted in decreasing (increasing) lexicographic
!    (y,x) order of coord; last NHOLA entries are for attached
!    holes; if bottom vertex of attached hole is a simple
!    vertex of boundary curve containing the hole then entry
!    contains index of bottom vertex otherwise entry contains
!    index of top vertex (which is simple).
!
!    Workspace, integer MAP(1:NCUR), used for mapping input boundary curve 
!    numbers to polygon numbers.
!
!    Workspace, HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG) - hash table and edge records
!    used to determine matching occurrences of separator or
!    cut interface edges by calling routine EDGHT.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 1, 2, 4, 5, 215, or 216.
!
  integer htsiz
  integer maxedg
  integer maxho
  integer maxhv
  integer maxpv
  integer ncur
!
  double precision angle
  integer edge(4,maxedg)
  integer, parameter :: edgv = 4
  logical first
  integer hdfree
  integer holv(maxho)
  integer ht(0:htsiz-1)
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer icur(ncur)
  integer ierror
  integer incr
  integer ipoly
  integer iv
  integer ivrt(*)
  integer ivs
  integer j
  integer j1
  integer j2
  integer jend
  integer jstr
  integer k
  integer kmax
  integer kmin
  integer kpoly
  integer l
  integer last
  integer, parameter :: loc = 1
  integer lv
  integer lvp
  integer lvs
  integer map(ncur)
  integer mpoly
  integer nh2
  integer nhola
  integer nhole
  integer nholi
  integer nht
  integer npolg
  integer nv
  integer nvbc(ncur)
  integer nvc
  integer nvert
  integer, parameter :: polg = 2
  integer pvl(4,maxpv)
  integer regnum(maxhv)
  integer, parameter :: succ = 3
  double precision vcl(2,nvc)
  double precision x
  double precision xmax
  double precision xmin
  double precision y
  double precision ymax
  double precision ymin
!
  ierror = 0
  nhola = 0
  nhole = 0
  nholi = 0
  nvert = 0

  do i = 1, ncur

    nvert = nvert + nvbc(i)

    if ( icur(i) > 0 ) then
      nhole = nhole + 1
    else if ( icur(i) < 0 ) then
      nholi = nholi + 1
      nvert = nvert + nvbc(i)
    end if

  end do

  npolg = ncur - nhole
  ipoly = 0
  iv = 0
  nv = 0
  hdfree = 0
  last = 0
  nht = 0

  ht(0:htsiz-1) = 0

  if ( ncur + nholi > maxhv ) then
    ierror = 4
    return
  else if ( nvert > maxpv ) then
    ierror = 5
    return
  else if ( ( nhole + nholi ) * 2 > maxho ) then
    ierror = 2
    return
  end if
!
!  Initialize REGNUM, HVL, PVL arrays for outer boundary curves.
!
  do i = 1, ncur

    if ( icur(i) /= 0 ) then
      map(i) = 0
      go to 40
    end if

    ipoly = ipoly + 1
    regnum(ipoly) = i
    hvl(ipoly) = iv + 1
    map(i) = ipoly
    jstr = nv + 1
    jend = nv + nvbc(i)

    do j = jstr, jend

      iv = iv + 1
      pvl(loc,iv) = abs ( ivrt(j) )
      pvl(polg,iv) = ipoly
      pvl(succ,iv) = iv + 1

      if ( ivrt(j) > 0 ) then
        pvl(edgv,iv) = 0
      else
!
!  The edge originating from current vertex is on a cut or
!  separator interface. Search in hash table for edge, and
!  insert or delete edge.  Set EDGV value if possible.
!
         lv = abs ( ivrt(j) )
         if ( lv > incr ) then
           lv = lv - incr
         end if

         if ( j < jend ) then
           lvs = abs ( ivrt(j+1) )
         else
           lvs = abs ( ivrt(jstr) )
         end if

         if ( lvs > incr ) then
           lvs = lvs - incr
         end if

         call edght ( lv, lvs, iv, nvc, htsiz, maxedg, hdfree, last, ht, &
           edge, ivs, ierror )

         if ( ierror /= 0 ) then
           return
         end if

         if ( ivs > 0 ) then
           pvl(edgv,iv) = ivs
           pvl(edgv,ivs) = iv
           nht = nht - 1
         else
           nht = nht + 1
         end if

       end if

     end do

     pvl(succ,iv) = hvl(ipoly)

40   continue

     nv = nv + nvbc(i)

  end do

  if ( nht /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSPGDC - Fatal error!'
    write ( *, '(a)' ) '  NHT /= 0.'
    ierror = 215
    return
  end if
!
!  Initialize REGNUM, HVL, PVL arrays for the hole interfaces.
!
  if ( nholi == 0 ) then
    go to 100
  end if

  do i = 1, ncur

    if ( icur(i) < 0 ) then
      ipoly = ipoly + 1
      map(i) = ipoly
    end if

  end do

  nv = 0

  do i = 1, ncur

    if ( icur(i) < 0 ) then

      ipoly = ipoly + 1
      kpoly = ipoly - nholi
      mpoly = map(-icur(i))
      regnum(kpoly) = i
      hvl(kpoly) = iv + 1
      hvl(ipoly) = iv + 2
      jstr = nv + 1
      jend = nv + nvbc(i)

      do j = jstr, jend
        iv = iv + 2
        pvl(loc,iv-1) = ivrt(j)
        pvl(polg,iv-1) = kpoly
        pvl(succ,iv-1) = iv + 1
        pvl(edgv,iv-1) = iv + 2
        pvl(loc,iv) = ivrt(j)
        pvl(polg,iv) = mpoly
        pvl(succ,iv) = iv - 2
        pvl(edgv,iv) = iv - 3
      end do

      pvl(succ,iv-1) = hvl(kpoly)
      pvl(edgv,iv-1) = hvl(ipoly)
      pvl(succ,hvl(ipoly)) = iv
      pvl(edgv,hvl(ipoly)) = iv - 1

    end if

    nv = nv + nvbc(i)
   
  end do
!
!  Initialize HVL, PVL arrays for the ordinary holes.
!
100 continue

  if ( nhole == 0 ) then
    go to 140
  end if

  nv = 0

  do i = 1, ncur

    if ( icur(i) > 0 ) then

      ipoly = ipoly + 1
      mpoly = map(icur(i))
      hvl(ipoly) = iv + 1
      jstr = nv + 1
      jend = nv + nvbc(i)

      do j = jstr, jend
        iv = iv + 1
        pvl(loc,iv) = ivrt(j)
        pvl(polg,iv) = mpoly
        pvl(succ,iv) = iv - 1
        pvl(edgv,iv) = 0
      end do

      pvl(succ,hvl(ipoly)) = iv

    end if

    nv = nv + nvbc(i)

  end do
!
!  Determine bottom or top simple vertex of attached holes.
!
  140 continue

  nhole = nhole + nholi
  nh2 = nhole + nhole
  j1 = 0
  j2 = 0

  do i = 1, npolg-nholi

   j = hvl(i)

  150    continue

      if ( pvl(loc,j) > incr ) then
         j = pvl(succ,j)
         if ( j /= hvl(i) ) then
          go to 150
         else
          ierror = 216
          return
         end if
      end if

   first = .true.

  160    continue

      lv = pvl(loc,j)

      if ( j1 > 0 ) then
        if ( lv <= incr ) then
          j2 = j
        else if ( lv - incr == lvs ) then
          j2 = j
        else
          pvl(loc,j) = lv - incr
        end if
      else if ( lv > incr ) then
        j1 = j
        lvs = lv - incr
        pvl(loc,j) = lvs
      end if

      if ( j2 > 0 ) then
!
!  (Part of) hole starts at vertex J1 and ends at J2.
!
         if ( lv <= incr .and. lv /= lvs ) go to 180
         k = j1

  170          continue

          if ( k == j1 ) then
            kmin = k
            kmax = k
            xmin = vcl(1,lvs)
            ymin = vcl(2,lvs)
            xmax = xmin
            ymax = ymin

          else

            l = pvl(loc,k)
            x = vcl(1,l)
            y = vcl(2,l)

            if ( y < ymin .or. y == ymin .and. x < xmin ) then
              kmin = k
              xmin = x
              ymin = y
            else if ( y > ymax .or. y == ymax .and. x > xmax ) then
              kmax = k
              xmax = x
              ymax = y
            end if

          end if

          k = pvl(succ,k)
         if ( k /= j2 ) then
           go to 170
         end if

         if ( kmin == j1 ) then
           kmin = kmax
         end if

         nhola = nhola + 1

         if ( nh2 + nhola > maxho ) then
          ierror = 2
          return
         end if

         holv(nh2+nhola) = kmin

180      continue

         j1 = 0
         j2 = 0

         if ( lv > incr ) then
           j1 = j
           pvl(loc,j) = lvs
         end if

      end if
      j = pvl(succ,j)

   if ( first ) then
     first = .false.
     jend = j
     go to 160
   else if ( j /= jend ) then
     go to 160
   end if

  end do
!
!  Initialize the IANG array.
!
  do i = 1, npolg+nhole

    j = hvl(i)
    lvp = pvl(loc,j)
    iv = pvl(succ,j)
    lv = pvl(loc,iv)

200   continue

      ivs = pvl(succ,iv)
      lvs = pvl(loc,ivs)
      iang(iv) = angle ( vcl(1,lvp), vcl(2,lvp), vcl(1,lv), vcl(2,lv), &
        vcl(1,lvs), vcl(2,lvs) )

      if ( iv /= j ) then
        lvp = lv
        iv = ivs
        lv = lvs
        go to 200
      end if

  end do
!
!  Initialize HOLV array.
!
  if ( nhole > 0 ) then
    call holvrt ( nhole, vcl, hvl(npolg+1), pvl, holv )
  end if

  return
end subroutine




subroutine dtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )
!
!*******************************************************************************
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2-D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!    On abnormal return, IERROR is set to 8, 224, or 225.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NPT, the number of vertices.
!
!    Input, integer MAXST, the maximum size available for the STACK array; 
!    should be about NPT to be safe, but MAX(10,2*LOG2(NPT)) is usually enough.
!
!    Input, double precision VCL(2,NPT), the coordinates of the vertices.
!
!    Input/output, integer IND(NPT), the indices in VCL of the vertices 
!    to be triangulated.  On output, IND has been permuted by the sort.
!
!    Output, integer NTRI, the number of triangles in the triangulation; 
!    NTRI is equal to 2*NPT - NB - 2, where NB is the number of boundary 
!    vertices.
!
!    Output, integer TIL(3,NTRI), the nodes that make up each triangle.
!    The elements are indices of VCL.  The vertices of the triangles are 
!    in counter clockwise order.
!
!    Output, integer TNBR(3,NTRI), the triangle neighbor list.
!    Positive elements are indices of TIL; negative elements are used for links
!    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Workspace, integer STACK(MAXST), used for a stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer IERROR, an error flag, nonzero if an error occurred.
!
  integer maxst
  integer npt
!
  double precision cmax
  integer e
  integer i
  integer ierror
  integer ind(npt)
  integer j
  integer k
  integer l
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer m
  integer m1
  integer m2
  integer, parameter :: msglvl = 0
  integer n
  integer ntri
  integer redg
  integer rtri
  integer stack(maxst)
  integer t
  integer til(3,npt*2)
  integer tnbr(3,npt*2)
  double precision tol
  integer top
  double precision vcl(2,npt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  ierror = 0
!
!  Sort the vertices.
!
  call dhpsrt ( 2, npt, 2, vcl, ind )
!
!  Ensure that no two consecutive points are too close.
!
  m1 = ind(1)

  do i = 2, npt

    m = m1
    m1 = ind(i)

    k = 0
    do j = 1, 2

      cmax = max ( abs ( vcl(j,m) ), abs ( vcl(j,m1) ) )

      if ( abs ( vcl(j,m) - vcl(j,m1) ) > tol * cmax .and. cmax > tol ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      ierror = 224
      return
    end if

  end do
!
!  Take the first two points, M1 and M2, and find a suitable non-collinear
!  third, M.  All points between M2 and M are very close to collinear
!  with M1 and M2.
!
  m1 = ind(1)
  m2 = ind(2)
  j = 3

  do

    if ( j > npt ) then
      ierror = 225
      return
    end if

    m = ind(j)
    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1
  
  end do

  ntri = j - 2
!
!  Depending on the orientation of M1, M2, and M, set up the initial
!  triangle data.
!
  if ( lr == -1 ) then

    til(1,1) = m1
    til(2,1) = m2
    til(3,1) = m
    tnbr(3,1) = -3

    do i = 2, ntri

      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m1
      til(2,i) = m2
      til(3,i) = m
      tnbr(1,i-1) = -3 * i
      tnbr(2,i-1) = i
      tnbr(3,i) = i - 1

    end do

    tnbr(1,ntri) = -3 * ntri - 1
    tnbr(2,ntri) = -5
    ledg = 2
    ltri = ntri

  else

    til(1,1) = m2
    til(2,1) = m1
    til(3,1) = m
    tnbr(1,1) = -4

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m2
      til(2,i) = m1
      til(3,i) = m
      tnbr(3,i-1) = i
      tnbr(1,i) = -3 * i - 3
      tnbr(2,i) = i - 1
    end do

    tnbr(3,ntri) = -3 * ntri
    tnbr(2,1) = -3 * ntri - 2
    ledg = 2
    ltri = 1

  end if

  if ( msglvl == 4 ) then

    m2 = ind(1)
    write ( *, 600 ) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)

    do i = 2, j-1
      m1 = m2
      m2 = ind(i)
      write ( *,600) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
      write ( *,600) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    end do

  end if
!
!  Insert vertices one at a time from outside the convex hull, determine
!  the visible boundary edges, and apply diagonal edge swaps until
!  the Delaunay triangulation of the vertices (so far) is obtained.
!
  top = 0

  do i = j+1, npt

    if ( msglvl == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, 600 ) i
    end if

    m = ind(i)
    m1 = til(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = til(ledg+1,ltri)
    else
      m2 = til(1,ltri)
    end if

    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0D+00 )

    if ( lr > 0 ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tnbr(ledg,ltri)
      rtri = l / 3
      redg = mod ( l, 3 ) + 1
    end if

    call vbedg ( vcl(1,m), vcl(2,m), vcl, til, tnbr, ltri, ledg, rtri, redg )

    n = ntri + 1
    l = -tnbr(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( top > maxst ) then
        ierror = 8
        return
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
      end if

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    if ( msglvl == 4 ) then
      write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
    end if

    tnbr(ledg,ltri) = -3 * n - 1
    tnbr(2,n) = -3 * ntri - 2
    tnbr(3,ntri) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, maxst, ltri, ledg, vcl, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      return
    end if

  end do

  if ( msglvl == 4 ) then
    write ( *, '(i7)' ) npt + 1
  end if

600 format (1x,i7,4f15.7)

  return
end subroutine




subroutine dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )
!
!*******************************************************************************
!
!! DTRIW2 constructs an incremental Delaunay triangulation in 2D.
!
!
!  Purpose: 
!
!    Construct Delaunay triangulation of 2-D vertices using
!    incremental approach and diagonal edge swaps. Vertices are
!    inserted one at a time in order given by IND array. The initial
!    triangles created due to a new vertex are obtained by a walk
!    through the triangulation until location of vertex is known.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NPT, the number of 2-D points (vertices).
!
!    Input, integer MAXST, the maximum size available for STACK array; should 
!    be about NPT to be safe, but MAX(10,2*LOG2(NPT)) usually enough.
!
!    Input, double precision VCL(1:2,1:*), the coordinates of 2-D vertices.
!
!    Input, integer IND(1:NPT), indices in VCL of vertices to be triangulated;
!    vertices are inserted in order given by this array.
!
!    Output, integer NTRI, the number of triangles in triangulation; equal to
!    2*NPT - NB - 2 where NB = number of boundary vertices.
!
!    Output, integer TIL(1:3,1:NTRI), the triangle incidence list; elements 
!    are indices of VCL; vertices of triangles are in counter clockwise order.
!
!    Output, integer TNBR(1:3,1:NTRI), the triangle neighbor list; positive 
!    elements are indices of TIL; negative elements are used for links
!    of counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Workspace, integer STACK(1:MAXST), used for stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer IERROR, error flag. For abnormal return,
!    IERROR is set to 8, 224, 225, or 226.
!
  integer maxst
  integer npt
!
  integer bedg
  integer btri
  double precision cmax
  integer e
  integer em1
  integer ep1
  integer ntri
  integer i
  integer i3
  integer ierror
  integer ind(npt)
  integer j
  integer l
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer m
  integer m1
  integer m2
  integer m3
  integer, parameter :: msglvl = 0
  integer n
  integer redg
  integer rtri
  integer stack(maxst)
  integer t
  integer til(3,npt*2)
  integer tnbr(3,npt*2)
  integer top
  double precision tol
  double precision vcl(2,*)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine the initial triangle.
!
  m1 = ind(1)
  m2 = ind(2)

  do j = 1, 2
    cmax = max ( abs ( vcl(j,m1) ), abs ( vcl(j,m2) ) )
    if ( abs ( vcl(j,m1) - vcl(j,m2) ) > tol * cmax .and. cmax > tol ) then
      go to 20
    end if
  end do

  ierror = 224
  return

20 continue

  i3 = 3

30 continue

  if ( i3 > npt ) then
    ierror = 225
    return
  end if

  m = ind(i3)
  lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
    vcl(2,m2), 0.0D+00 )

  if ( lr == 0 ) then
    i3 = i3 + 1
    go to 30
  end if

  if ( i3 /= 3 ) then
    ind(i3) = ind(3)
    ind(3) = m
  end if

  ntri = 1

  if ( lr == -1 ) then
    til(1,1) = m1
    til(2,1) = m2
  else
    til(1,1) = m2
    til(2,1) = m1
  end if

  til(3,1) = m
  tnbr(1,1) = -4
  tnbr(2,1) = -5
  tnbr(3,1) = -3

  if ( msglvl == 4 ) then
    write ( *,600) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
    write ( *,600) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
  end if
!
!  Insert vertices one at a time from anywhere.
!  Walk through the triangulation to determine the location of the new vertex.
!  Apply diagonal edge swaps until Delaunay triangulation of vertices
!  (so far) is obtained.
!
  top = 0

  do i = 4, npt

    if ( msglvl == 4 ) then
      write ( *,600) i
    end if

    m = ind(i)
    rtri = ntri

    call walkt2 ( vcl(1,m), vcl(2,m), ntri, vcl, til, tnbr, rtri, redg, ierror )

    if ( redg == 0 ) then

      m1 = til(1,rtri)
      m2 = til(2,rtri)
      m3 = til(3,rtri)
      til(3,rtri) = m

      if ( tnbr(1,rtri) > 0 ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m2
      til(2,ntri) = m3
      til(3,ntri) = m
      n = tnbr(2,rtri)
      tnbr(1,ntri) = n

      if ( n > 0 ) then
        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if
        top = top + 1
        stack(top) = ntri
      end if

      ntri = ntri + 1
      til(1,ntri) = m3
      til(2,ntri) = m1
      til(3,ntri) = m
      n = tnbr(3,rtri)
      tnbr(1,ntri) = n

      if ( n > 0 ) then
        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if
        top = top + 1
        stack(top) = ntri
      end if

      tnbr(2,rtri) = ntri - 1
      tnbr(3,rtri) = ntri
      tnbr(2,ntri-1) = ntri
      tnbr(3,ntri-1) = rtri
      tnbr(2,ntri) = rtri
      tnbr(3,ntri) = ntri - 1

      if ( tnbr(1,ntri-1) <= 0 ) then

        t = rtri
        e = 1

        do

          if ( tnbr(e,t) <= 0 ) then
            exit
          end if

          t = tnbr(e,t)

          if ( til(1,t) == m2 ) then
            e = 3
          else if ( til(2,t) == m2 ) then
            e = 1
          else
            e = 2
          end if

        end do

        tnbr(e,t) = -3 * ntri + 3

      end if

      if ( tnbr(1,ntri) <= 0 ) then

        t = ntri - 1
        e = 1

        do

          if ( tnbr(e,t) <= 0 ) then
            exit
          end if

          t = tnbr(e,t)
          if ( til(1,t) == m3 ) then
            e = 3
          else if ( til(2,t) == m3 ) then
            e = 1
          else
            e = 2
          end if

        end do

        tnbr(e,t) = -3 * ntri

      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

   else if ( redg < 0 ) then

      redg = -redg
      ltri = 0
      call vbedg ( vcl(1,m), vcl(2,m), vcl, til, tnbr, ltri, ledg, rtri, redg )
      n = ntri + 1
      l = -tnbr(ledg,ltri)

60    continue

         t = l / 3
         e = mod ( l, 3 ) + 1
         l = -tnbr(e,t)
         m2 = til(e,t)

         if ( e <= 2 ) then
           m1 = til(e+1,t)
         else
           m1 = til(1,t)
         end if

         ntri = ntri + 1
         tnbr(e,t) = ntri
         til(1,ntri) = m1
         til(2,ntri) = m2
         til(3,ntri) = m
         tnbr(1,ntri) = t
         tnbr(2,ntri) = ntri - 1
         tnbr(3,ntri) = ntri + 1
         top = top + 1

         if ( top > maxst ) then
           ierror = 8
           go to 100
         end if
 
        stack(top) = ntri

         if ( msglvl == 4 ) then
           write (*,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
         end if

      if ( t /= rtri .or. e /= redg ) then
        go to 60
      end if

      if ( msglvl == 4 ) then
        write (*,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
      end if

      tnbr(ledg,ltri) = -3*n - 1
      tnbr(2,n) = -3*ntri - 2
      tnbr(3,ntri) = -l

   else if ( redg <= 3 ) then

      m1 = til(redg,rtri)

      if ( redg == 1 ) then
        e = 2
        ep1 = 3
      else if ( redg == 2 ) then
        e = 3
        ep1 = 1
      else
        e = 1
        ep1 = 2
      end if

      m2 = til(e,rtri)
      til(e,rtri) = m
      m3 = til(ep1,rtri)

      if ( tnbr(ep1,rtri) > 0 ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m
      til(2,ntri) = m2
      til(3,ntri) = m3
      n = tnbr(e,rtri)
      tnbr(2,ntri) = n
      tnbr(3,ntri) = rtri
      tnbr(e,rtri) = ntri

      if ( n > 0 ) then
        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if
        top = top + 1
        stack(top) = ntri
      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

      ltri = tnbr(redg,rtri)

      if ( ltri <= 0 ) then
        tnbr(1,ntri) = ltri
        tnbr(redg,rtri) = -3*ntri
        if ( tnbr(2,ntri) <= 0 ) then
          tnbr(1,ntri) = -3*ntri - 1
        end if

      else

        tnbr(1,ntri) = ntri + 1
        tnbr(redg,rtri) = ltri

        if ( til(1,ltri) == m2 ) then
          ledg = 1
          em1 = 2
          e = 3
        else if ( til(2,ltri) == m2 ) then
          ledg = 2
          em1 = 3
          e = 1
        else
          ledg = 3
          em1 = 1
          e = 2
        end if

        til(ledg,ltri) = m
        m3 = til(e,ltri)
        if ( tnbr(em1,ltri) > 0 ) then
          top = top + 1
          stack(top) = ltri
        end if
        ntri = ntri + 1
        til(1,ntri) = m2
        til(2,ntri) = m
        til(3,ntri) = m3
        tnbr(1,ntri) = ntri - 1
        tnbr(2,ntri) = ltri
        n = tnbr(e,ltri)
        tnbr(3,ntri) = n
        tnbr(e,ltri) = ntri

        if ( n > 0 ) then
          if ( tnbr(1,n) == ltri ) then
            tnbr(1,n) = ntri
          else if ( tnbr(2,n) == ltri ) then
            tnbr(2,n) = ntri
          else
            tnbr(3,n) = ntri
          end if
          top = top + 1
          stack(top) = ntri
        end if

        if ( msglvl == 4 ) then
          write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
        end if

        if ( tnbr(2,ntri-1) <= 0 ) then

          t = ntri
          e = 3

          do

            if ( tnbr(e,t) <= 0 ) then
              exit
            end if

            t = tnbr(e,t)
            if ( til(1,t) == m2 ) then
              e = 3
            else if ( til(2,t) == m2 ) then
              e = 1
            else
              e = 2
            end if

          end do

          tnbr(e,t) = -3 * ntri + 2

        end if

         if ( tnbr(3,ntri) <= 0 ) then

          t = ltri

          if ( ledg <= 2 ) then
             e = ledg + 1
          else
             e = 1
          end if

          do

            if ( tnbr(e,t) <= 0 ) then
              exit
            end if

            t = tnbr(e,t)
            if ( til(1,t) == m3 ) then
              e = 3
            else if ( til(2,t) == m3 ) then
              e = 1
            else
              e = 2
            end if

          end do

          tnbr(e,t) = -3 * ntri - 2

         end if

      end if

    else
      ierror = 224
      go to 100
    end if

    btri = 0
    bedg = 0

    call swapec ( m, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

  end do

100 continue

  if ( i3 /= 3 ) then
    t = ind(i3)
    ind(i3) = ind(3)
    ind(3) = t
  end if

  if ( msglvl == 4 ) then
    write ( *,600) npt+1
  end if

600 format (1x,i7,4f15.7)

  return
end subroutine



subroutine edght ( a, b, v, n, htsiz, maxedg, hdfree, last, ht, edge, w, ierror )
!
!*******************************************************************************
!
!! EDGHT searches a hash table for a record in EDGE containing key (A,B).
!
!
!  Purpose: 
!
!    Search in hash table HT for record in EDGE containing
!    key (A,B).
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer A, B, vertex indices, > 0, of edge (also key of hash table).
!
!    Input, integer V, value associated with edge.
!
!    Input, integer N, upper bound on A, B.
!
!    Input, integer HTSIZ, the size of hash table HT.
!
!    Input, integer MAXEDG, the maximum size available for EDGE array.
!
!    Input/output, integer HDFREE, head pointer to linked list of free entries 
!    of EDGE array due to deletions.  Before first call to this routine, HDFREE
!    should be set to 0.
!
!    Input/output, integer LAST, index of last entry used in EDGE array.  Before 
!    first call to this routine, LAST should be set to 0.
!
!    Input/output, integer HT(0:HTSIZ-1), hash table of head pointers (direct 
!    chaining with ordered lists is used).  Before first call to this routine, 
!    entries of HT should be set to 0.  If key with A,B is found then this 
!    record is deleted from hash table, else record is inserted in hash table.
!
!    Input/output, integer EDGE(1:4,1:MAXEDG), entries of hash table records;
!      EDGE(1,I) = MIN(A,B); EDGE(2,I) = MAX(A,B);
!      EDGE(3,I) = V; EDGE(4,I) = link
!    If key with A,B is found then this record is deleted
!    from hash table, else record is inserted in hash table.
!
!    Output, integer W, EDGE(3,INDEX), where INDEX is index of record, if found;
!    else 0.
!
!    Output, integer IERROR, error flag.  For abnormal return,
!    IERROR is set to 1
!
  integer htsiz
  integer maxedg
!
  integer a
  integer aa
  integer b
  integer bb
  integer bptr
  integer edge(4,maxedg)
  integer hdfree
  integer ht(0:htsiz-1)
  integer ierror
  integer k
  integer last
  integer n
  integer newp
  integer ptr
  integer v
  integer w
!
  ierror = 0

  if ( a < b ) then
    aa = a
    bb = b
  else
    aa = b
    bb = a
  end if

  k = mod ( aa*n + bb, htsiz )
  bptr = -1
  ptr = ht(k)

10 continue

  if ( ptr /= 0 ) then

    if ( edge(1,ptr) > aa ) then

      go to 20

    else if ( edge(1,ptr) == aa ) then

      if ( edge(2,ptr) > bb ) then

        go to 20

      else if ( edge(2,ptr) == bb ) then

        if ( bptr == -1 ) then
          ht(k) = edge(4,ptr)
        else
          edge(4,bptr) = edge(4,ptr)
        end if

        edge(4,ptr) = hdfree
        hdfree = ptr
        w = edge(3,ptr)
        return

      end if

    end if

    bptr = ptr
    ptr = edge(4,ptr)
    go to 10

  end if

20 continue

  if ( hdfree > 0 ) then

    newp = hdfree
    hdfree = edge(4,hdfree)

  else

    last = last + 1
    newp = last

    if ( last > maxedg ) then
      ierror = 1
      return
    end if

  end if

  if ( bptr == -1 ) then
    ht(k) = newp
  else
    edge(4,bptr) = newp
  end if

  edge(1,newp) = aa
  edge(2,newp) = bb
  edge(3,newp) = v
  edge(4,newp) = ptr
  w = 0

  return
end subroutine


subroutine eqdis2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid,  &
  nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl,  &
  pvl, iang, area, psi, h, iwk, wk, ierror )
!
!*******************************************************************************
!
!! EQDIS2 further subdivides convex polygons for mesh equidistribution.
!
!
!  Purpose: 
!
!    Further subdivide convex polygons so that an approximately
!    equidistributing triangular mesh can be constructed with
!    respect to a heuristic or a user-supplied mesh distribution
!    function (MDF), and determine triangle size for each polygon of
!    decomposition.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical HFLAG, is .TRUE. if heuristic mdf, .FALSE. if
!    user-supplied mdf.
!
!    Input, external UMDF, a user-supplied mesh distribution function 
!    of the form:
!
!      function umdf ( x, y )
!      double precision umdf
!      double precision x
!      double precision y
!
!    Input, double precision KAPPA, the mesh smoothness parameter in 
!    the interval [0.0,1.0].
!
!    Input, double precision ANGSPC, the angle spacing parameter in radians
!    used to determine extra points as possible endpoints of separators.
!
!    Input, double precision ANGTOL, the angle tolerance parameter in radians
!    used in accepting separators.
!
!    Input, double precision DMIN, a parameter used to determine if variation 
!    of mdf in polygon is 'sufficiently high'.
!
!    Input, integer NMIN, a parameter used to determine if 'sufficiently large'
!    number of triangles in polygon.
!
!    Input, integer NTRID, the desired number of triangles in mesh.
!
!    Input/output, integer NVC, the number of vertex coordinates or 
!    positions used in VCL array.
!
!    Input/output, integer NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array, should 
!    be >= number of vertex coordinates required for decomposition
!    (approximately NVC + 2*NS where NS is expected number of new separators).
!
!    Input, integer MAXHV, the maximum size available for HVL, REGNUM, AREA, 
!    PSI, H arrays; should be >= number of polygons required for
!    decomposition (approximately NPOLG + NS).
!
!    Input, integer MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be >= number of polygon vertices required for decomposition
!    (approximately NVERT + 5*NS).
!
!    Input, integer MAXIW, the maximum size available for IWK array; should 
!    be >= MAX(2*NP, NVERT + NPOLG + 3*NVRT + INT(2*PI/ANGSPC))
!    where NVRT is maximum number of vertices in a convex
!    polygon of the (input) decomposition, NP is expected
!    value of NPOLG on output.
!
!    Input, integer MAXWK, the maximum size available for WK array; should
!    be >= NVC + NVERT + 2*NPOLG + 3*(NVRT + INT(2*PI/ANGSPC)).
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, double precision HVL(1:NPOLG), head vertex list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT), the 
!    polygon vertex list and interior angles; see routine DSPGDC for more 
!    details.  Note: The data structures should be as output from routine
!    CVDEC2.
!
!    Output, double precision AREA(1:NPOLG), the area of convex polygons 
!    in the decomposition.
!
!    Output, double precision PSI(1:NPOLG), the smoothed mean mdf 
!    values in the convex polygons.
!
!    Output, double precision H(1:NPOLG), the triangle size for convex polygons.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  For abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 200, 201, or 222.
!
  integer maxhv
  integer maxiw
  integer maxpv
  integer maxvc
  integer maxwk
!
  double precision angspc
  double precision angtol
  double precision area(maxhv)
  double precision dmin
  integer edgval
  double precision h(maxhv)
  logical hflag
  integer hvl(maxhv)
  double precision iang(maxpv)
  integer ierror
  integer ivrt
  integer iwk(maxiw)
  double precision kappa
  integer m
  integer n
  integer nmin
  integer npolg
  integer ntrid
  integer nvc
  integer nvert
  double precision psi(maxhv)
  integer pvl(4,maxpv)
  integer regnum(maxhv)
  double precision umdf
  double precision vcl(2,maxvc)
  integer vrtval
  integer widsq
  double precision wk(maxwk)
  integer xivrt
!
  external umdf
!
  ierror = 0
  ivrt = 1
  xivrt = ivrt + nvert
  m = xivrt + npolg

  if ( m > maxiw ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  M > MAXIW'
    ierror = 6
    return
  end if

  widsq = 1

  if ( hflag ) then

    edgval = widsq + npolg
    vrtval = edgval + nvert
    n = npolg + nvert + nvc

    if ( n > maxwk ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
      write ( *, '(a)' ) '  N > MAXWK'
      ierror = 7
      return
    end if

  else

    edgval = 1
    vrtval = 1
    n = 0

  end if

  call dsmdf2 ( hflag, nvc, npolg, maxwk-n, vcl, hvl, pvl, iang, iwk(ivrt), &
    iwk(xivrt), wk(widsq), wk(edgval), wk(vrtval), area, wk(n+1), ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  DSMDF2 returns error condition.'
    return
  end if

  call mfdec2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid, nvc, &
    npolg, nvert, maxvc, maxhv, maxpv, maxiw-m, maxwk-n, vcl, regnum, hvl, &
    pvl, iang, iwk(ivrt), iwk(xivrt), wk(widsq), wk(edgval), wk(vrtval), &
    area, psi, iwk(m+1), wk(n+1), ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  MFDEC2 returns error condition.'
    return
  end if

  if ( 2 * npolg > maxiw ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EQDIS2 - Fatal error!'
    write ( *, '(a)' ) '  2 * NPOLG > MAXIW.'
    ierror = 6
    return
  end if

  call trisiz ( ntrid, npolg, hvl, pvl, area, psi, h, iwk, iwk(npolg+1) )

  return
end subroutine



subroutine fndsep ( angac1, xr, yr, nvrt, xc, yc, ivis, theta, nv, iv,  &
  vcl, pvl, iang, angsep, i1, i2, wkang )
!
!*******************************************************************************
!
!! FNDSEP finds separators to resolve a reflex vertex.
!
!
!  Purpose: 
!
!    Find 1 or 2 separators which can resolve a reflex vertex
!    (XR,YR) using a max-min angle criterion from list of vertices
!    in increasing polar angle with respect to the reflex vertex.
!
!    Preference is given to 1 separator.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGAC1, the angle tolerance parameter used 
!    for preference in accepting one separator.
!
!    Input, double precision XR, YR, the coordinates of reflex vertex.
!
!    Input, integer NVRT, (number of vertices) - 1.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    of possible endpoints of a separator.
!
!    Input, integer IVIS(0:NVRT), contains information about the vertices of
!    XC, YC arrays with respect to the polygon vertex list; if
!    IVIS(I) > 0 then vertex (XC(I),YC(I)) has index IVIS(I)
!    in PVL; if IVIS(I) < 0 then vertex (XC(I),YC(I)) is on
!    the edge joining vertices with indices -IVIS(I) and
!    SUCC(-IVIS(I)) in PVL.
!
!    Input, double precision THETA(0:NVRT), the polar angles of vertices 
!    in increasing order; THETA(NVRT) is the interior angle of reflex vertex;
!    THETA(I), I >= 0, is the polar angle of (XC(I),YC(I))
!    with respect to reflex vertex.
!
!    Input, integer NV, (number of vertices to be considered as endpoint of a
!    separator) - 1.
!
!    Input, integer IV(0:NV), the indices of vertices in XC, YC arrays to be
!    considered as endpoint of a separator; angle between
!    consecutive vertices is assumed to be < 180 degrees.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer PVL(1:4,1:*), double precision IANG(1:*), the polygon 
!    vertex list, interior angles.
!
!    Output, double precision ANGSEP, the minimum of the 4 or 7 angles at the
!    boundary resulting from 1 or 2 separators, respectively.
!
!    Output, integer I1, I2, the indices of endpoints of separators in XC, 
!    YC arrays; I2 = -1 if there is only one separator, else I1 < I2.
!
!    Workspace, double precision WKANG(0:NV).
!
  double precision ang
  double precision angac1
  double precision angsep
  double precision angsp2
  integer i
  integer i1
  integer i2
  double precision iang(*)
  integer ii
  integer k
  integer l
  integer m
  integer nl
  integer nr
  integer nv
  integer nvrt
  integer iv(0:nv)
  integer ivis(0:nvrt)
  double precision minang
  integer p
  double precision phi
  double precision d_pi
  integer pvl(4,*)
  integer q
  integer r
  double precision theta(0:nvrt)
  double precision tol
  double precision vcl(2,*)
  double precision wkang(0:nv)
  double precision xc(0:nvrt)
  double precision xr
  double precision yc(0:nvrt)
  double precision yr
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine the vertices in the inner cone - indices P to Q.
!
  i = 0
  p = -1
  phi = theta(nvrt) - d_pi ( ) + tol

  do while ( p < 0 )

    if ( theta(iv(i)) >= phi ) then
      p = i
    else
      i = i + 1
    end if

  end do

  i = nv
  q = -1
  phi = d_pi ( ) - tol

  do while ( q < 0 )

    if ( theta(iv(i)) <= phi ) then
      q = i
    else
      i = i - 1
    end if

  end do
!
!  Use the max-min angle criterion to find the best separator
!  in inner cone.
!
  angsep = 0.0

  do i = p, q

    k = iv(i)
    ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
      vcl, pvl, iang )

    if ( ang > angsep ) then
      angsep = ang
      ii = iv(i)
    end if

  end do

  angsp2 = angsep
  if ( angsep >= angac1 ) then
    go to 110
  end if
!
!  If the best separator in inner cone is not 'good' enough,
!  use max-min angle criterion to try to find a better pair
!  of separators from the right and left cones.
!
  nr = 0
  nl = 0

  do r = 0, p-1

    wkang(r) = 0.0D+00

    if ( theta(iv(r)) > angsep ) then

      k = iv(r)

      ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
        vcl, pvl, iang )

      if ( ang > angsep ) then
        nr = nr + 1
        wkang(r) = ang
      end if

    end if

  end do

  if ( nr == 0 ) then
    go to 110
  end if

  phi = theta(nvrt) - angsep

  do l = q+1, nv

    wkang(l) = 0.0D+00

    if ( theta(iv(l)) < phi ) then

      k = iv(l)
      ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
        vcl, pvl, iang )

      if ( ang > angsep ) then
        nl = nl + 1
        wkang(l) = ang
      end if

    end if

  end do

  if ( nl == 0 ) then
    go to 110
  end if
!
!  Check all possible pairs for the best pair of separators
!  in the right and left cones.
!
  m = nv

  do r = p-1, 0, -1

     if ( m > q .and. wkang(r) > angsp2 ) then

      phi = theta(iv(r))

80    continue

      if ( m > q .and. ( wkang(m) <= angsp2 .or. &
        theta(iv(m)) - phi > d_pi ( ) - tol) ) then
        m = m - 1
        go to 80
      end if

      do l = q+1, m

         if ( wkang(l) > angsp2 ) then

            ang = min ( theta(iv(l)) - phi, wkang(r), wkang(l) )

            if ( ang > angsp2 ) then
             angsp2 = ang
             i1 = iv(r)
             i2 = iv(l)
            end if

         end if

       end do

     end if

  end do
!
!  Choose 1 or 2 separators based on max-min angle criterion or
!  ANGAC1 parameter.
!
  110 continue

  if ( angsp2 <= angsep ) then
    i1 = ii
    i2 = -1
  else
    angsep = angsp2
  end if

  return
end subroutine



subroutine fndtri ( iedg, mxtr, sflag, tedg, itr, ind, ierror )
!
!*******************************************************************************
!
!! FNDTRI finds two triangles containing a given edge.
!
!
!  Purpose: 
!
!    Find two triangles containing edge with index IEDG in array TEDG.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer IEDG, the index of edge to be searched in TEDG.
!
!    Input, integer MXTR, the maximum index of triangle to be searched in TEDG.
!
!    Input, logical SFLAG, is .TRUE. if and only if the second triangle is to be 
!    searched from end of array.
!
!    Input, integer TEDG(1:3,1:MXTR), triangle edge indices; see routine CVDTRI.
!
!    Output, integer ITR(1:2), IND(1:2), indices such that IEDG =
!    TEDG(IND(1),ITR(1)) = TEDG(IND(2),ITR(2)).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 231.
!
  integer mxtr
!
  integer i
  integer iedg
  integer ierror
  integer ind(2)
  integer itr(2)
  integer j
  integer k
  logical sflag
  integer tedg(3,mxtr)
!
!  Search from end of array TEDG.
!
  ierror = 0
  k = 1
  j = 1
  i = mxtr

10 continue

  do

    if ( tedg(j,i) == iedg ) then
      exit
    end if

    j = j + 1

    if ( j > 3 ) then
      j = 1
      i = i - 1
      if ( i <= 0 ) then
        ierror = 231
        return
      end if
    end if

  end do

  itr(k) = i
  ind(k) = j

  if ( k == 2 ) then
    return
  end if

  k = 2

  if ( sflag ) then

    j = 1
    i = i - 1

    if ( i <= 0 ) then
      ierror = 231
      return
    end if

    go to 10

  end if
!
!  Search from beginning of array TEDG for second triangle.
!
  j = 1
  i = 1

20 continue

  if ( i >= itr(1) ) then
    ierror = 231
    return
  end if

30 continue

  if ( tedg(j,i) /= iedg ) then
    j = j + 1
    if ( j > 3 ) then
      j = 1
      i = i + 1
      go to 20
    else
      go to 30
    end if
  end if

  itr(2) = i
  ind(2) = j

  return
end
subroutine gtime ( time )
!
!*******************************************************************************
!
!! GTIME gets the current CPU time in seconds.
!
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real TIME, the current reading of the CPU clock in seconds.
!
  integer clock_count
  integer clock_max
  integer clock_rate
  real time
!
  call system_clock ( clock_count, clock_rate, clock_max )

  time = real ( clock_count ) / real ( clock_rate )

  return
end subroutine



subroutine holvrt ( nhole, vcl, hvl, pvl, holv )
!
!*******************************************************************************
!
!! HOLVRT determines top and bottom vertices of holes in polygonal regions.
!
!
!  Purpose: 
!
!    Determine top and bottom vertices of holes in polygonal
!    regions, and sort top vertices in decreasing (y,x) order
!    and bottom vertices in increasing (y,x) order.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NHOLE, the number of holes in region(s).
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer HVL(1:NHOLE), the head vertex list; HVL(I) is index in PVL of
!    head vertex of Ith hole.
!
!    Input, integer PVL(1:4,1:*), the polygon vertex list; see routine DSPGDC.
!
!    Output, integer HOLV(1:NHOLE*2), the indices in PVL of top and bottom 
!    vertices of holes; first (last) NHOLE entries are for top (bottom)
!    vertices; top (bottom) vertices are sorted in decreasing
!    (increasing) lexicographic (y,x) order of coordinates.
!
  integer nhole
!
  integer holv(nhole*2)
  integer hv
  integer hvl(nhole)
  integer i
  integer imax
  integer imin
  integer iv
  integer j
  integer, parameter :: loc = 1
  integer lv
  integer nhp1
  integer pvl(4,*)
  integer, parameter :: succ = 3
  double precision vcl(2,*)
  double precision x
  double precision xmax
  double precision xmin
  double precision y
  double precision ymax
  double precision ymin
!
!  Determine top and bottom vertices of holes.
!
  do i = 1, nhole

    hv = hvl(i)
    iv = hv

    do

      lv = pvl(loc,iv)

      if ( iv == hv ) then

        imin = iv
        imax = iv
        xmin = vcl(1,lv)
        ymin = vcl(2,lv)
        xmax = xmin
        ymax = ymin

      else

        x = vcl(1,lv)
        y = vcl(2,lv)

        if ( y < ymin .or. y == ymin .and. x < xmin ) then
          imin = iv
          xmin = x
          ymin = y
        else if ( y > ymax .or. y == ymax .and. x  >  xmax ) then
          imax = iv
          xmax = x
          ymax = y
        end if

      end if

      iv = pvl(succ,iv)

      if ( iv == hv) then
        exit
      end if

    end do

    holv(i) = imax
    holv(i+nhole) = imin

  end do
!
!  Use linear insertion sort to sort the top vertices of holes
!  in decreasing (y,x) order, then bottom vertices in increasing
!  (y,x) order. It is assumed NHOLE is small.
!
  do i = 2, nhole

    hv = holv(i)
    lv = pvl(loc,hv)
    x = vcl(1,lv)
    y = vcl(2,lv)
    j = i

30  continue

    iv = holv(j-1)
    lv = pvl(loc,iv)

    if ( y > vcl(2,lv) .or. y == vcl(2,lv) .and. x > vcl(1,lv) ) then
      holv(j) = iv
      j = j - 1
      if ( j > 1) go to 30
    end if

    holv(j) = hv

  end do

  nhp1 = nhole + 1

  do i = nhp1+1, nhole+nhole

    hv = holv(i)
    lv = pvl(loc,hv)
    x = vcl(1,lv)
    y = vcl(2,lv)
    j = i

50  continue

    iv = holv(j-1)
    lv = pvl(loc,iv)

    if ( y < vcl(2,lv) .or. y == vcl(2,lv) .and. x < vcl(1,lv) ) then
      holv(j) = iv
      j = j - 1
      if ( j > nhp1) go to 50
    end if

    holv(j) = hv

  end do

  return
end subroutine


function i_modp ( i, j )
!
!*******************************************************************************
!
!! I_MODP returns the nonnegative remainder of integer division.
!
!
!  Formula:
!
!    If
!      NREM = I_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I_MODP(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I     J     MOD  I_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I_MODP, the nonnegative remainder when I is
!    divided by J.
!
  integer i
  integer i_modp
  integer j
!
  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I_MODP ( I, J ) called with J = ', j
    stop
  end if

  i_modp = mod ( i, j )

  if ( i_modp < 0 ) then
    i_modp = i_modp + abs ( j )
  end if

  return
end function



subroutine i_swap ( i, j )
!
!*******************************************************************************
!
!! I_SWAP swaps two integer values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  integer i
  integer j
  integer k
!
  k = i
  i = j
  j = k

  return
end subroutine


function i_wrap ( ival, ilo, ihi )
!
!*******************************************************************************
!
!! I_WRAP forces an integer to lie between given limits by wrapping.
!
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer I_WRAP, a "wrapped" version of IVAL.
!
  integer i_modp
  integer i_wrap
  integer ihi
  integer ilo
  integer ival
  integer wide
!
  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i_wrap = ilo
  else
    i_wrap = ilo + i_modp ( ival-ilo, wide )
  end if

  return
end function


subroutine ihpsrt ( k, n, lda, a, map )
!
!*******************************************************************************
!
!! IHPSRT uses heapsort on integer points in K-dimension.
!
!
!  Purpose:
!
!    Use heapsort to obtain the permutation of N K-dimensional
!    integer points so that the points are in lexicographic
!    increasing order.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer K, the dimension of points.
!
!    Input, integer N, the number of points.
!
!    Input, integer LDA, the leading dimension of array A in calling 
!    routine; should be >= K.
!
!    Input, integer A(1:K,1:*), the array of >= N K-dimensional integer points.
!
!    Input/output, integer MAP(1:N), the points of A with indices MAP(1), 
!    MAP(2), ..., MAP(N) are to be sorted.
!    On output, elements are permuted so that A(*,MAP(1)) <=
!    A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  integer lda
  integer n
!
  integer a(lda,*)
  integer i
  integer k
  integer map(n)
  integer t
!
  do i = n/2, 1, -1
    call isftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call isftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end subroutine

 
function iless ( k, p, q )
!
!*******************************************************************************
!
!! ILESS determines whether a K-dimensional point P is lexically less than Q.
!
!
!  Purpose: 
!
!    Determine whether P is lexicographically less than Q.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer K, the dimension of points.
!
!    Input, integer P(1:K), Q(1:K), two points to compare.
!
!    Output, logical ILESS, is .TRUE. if P < Q, .FALSE. otherwise.
!
  integer k
!
  integer i
  logical iless
  integer p(k)
  integer q(k)
!
  do i = 1, k

    if ( p(i) /= q(i) ) then

      if ( p(i) < q(i) ) then
        iless = .true.
      else
        iless = .false.
      end if

      return

    end if

  end do

  iless = .false.

  return
end function



subroutine imat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! IMAT_PRINT prints an integer matrix.
!
!
!  Modified:
!
!    08 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, integer A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer j
  integer jhi
  integer jlo
  integer m
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end subroutine


subroutine insed2 ( v, w, npolg, nvert, maxhv, maxpv, vcl, regnum, &
  hvl, pvl, iang, ierror )
!
!*******************************************************************************
!
!! INSED2 inserts an edge into the head and polygon vertex lists.
!
!
!  Purpose: 
!
!    Insert edge joining vertices V, W into head vertex
!    list and polygon vertex list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer V, W, indices in PVL of vertices which are the endpoints
!    of an edge to be added to decomposition.
!
!    Input, integer NPOLG, the number of positions used in HVL array.
!
!    Input, integer NVERT, the number of positions used in PVL array.
!
!    Input, integer MAXHV, the maximum size available for HVL array.
!
!    Input, integer MAXPV, the maximum size available for PVL array.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT), 
!    the polygon vertex list and interior angles.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 4 or 5.
!
  integer maxhv
  integer maxpv
!
  double precision angle
  integer, parameter :: edgv = 4
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer ierror
  integer l
  integer, parameter :: loc = 1
  integer lv
  integer lw
  integer, parameter :: msglvl = 0
  integer npolg
  integer nvert
  integer, parameter :: polg = 2
  integer pvl(4,maxpv)
  integer regnum(maxhv)
  integer, parameter :: succ = 3
  integer v
  double precision vcl(2,*)
  integer vv
  integer w
  integer ww
!
  ierror = 0

  if ( npolg >= maxhv ) then
    ierror = 4
    return
  else if ( nvert+2 > maxpv ) then
    ierror = 5
    return
  end if
!
!  Split linked list of vertices of the polygon containing vertices
!  V and W into two linked list of vertices of polygons with common
!  edge joining V and W.
!
  nvert = nvert + 2
  vv = nvert - 1
  ww = nvert
  lv = pvl(loc,v)
  lw = pvl(loc,w)
  pvl(loc,vv) = lv
  pvl(loc,ww) = lw
  pvl(polg,ww) = pvl(polg,v)
  pvl(succ,vv) = pvl(succ,v)
  pvl(succ,ww) = pvl(succ,w)
  pvl(succ,v) = ww
  pvl(succ,w) = vv
  pvl(edgv,vv) = pvl(edgv,v)
  pvl(edgv,ww) = pvl(edgv,w)
  pvl(edgv,v) = w
  pvl(edgv,w) = v

  if ( pvl(edgv,vv) > 0 ) then
    pvl(edgv,pvl(edgv,vv)) = vv
  end if

  if ( pvl(edgv,ww) > 0 ) then
    pvl(edgv,pvl(edgv,ww)) = ww
  end if

  l = pvl(loc,pvl(succ,vv))

  iang(vv) = angle ( vcl(1,lw), vcl(2,lw), vcl(1,lv), vcl(2,lv), vcl(1,l), &
    vcl(2,l) )

  iang(v) = iang(v) - iang(vv)
  l = pvl(loc,pvl(succ,ww))
  iang(ww) = angle ( vcl(1,lv), vcl(2,lv), vcl(1,lw), vcl(2,lw), vcl(1,l), &
    vcl(2,l) )

  iang(w) = iang(w) - iang(ww)
  npolg = npolg + 1
  i = vv

  do

    pvl(polg,i) = npolg
    i = pvl(succ,i)

    if ( i == vv ) then
      exit
    end if

  end do

  hvl(pvl(polg,v)) = v
  hvl(npolg) = vv
  regnum(npolg) = regnum(pvl(polg,v))

  if ( msglvl == 2 ) then
    write ( *,600) v,w,vcl(1,lv),vcl(2,lv),vcl(1,lw),vcl(2,lw)
  end if

600 format (1x,2i7,4f15.7)

  return
end subroutine



subroutine insvr2 ( xi, yi, wp, nvc, nvert, maxvc, maxpv, vcl, pvl, &
  iang, w, ierror )
!
!*******************************************************************************
!
!! INSVR2 inserts a point into the vertex coordinate and polygon vertex lists.
!
!
!  Purpose: 
!
!    Insert point (XI,YI) into vertex coordinate list and
!    polygon vertex list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XI, YI, the coordinates of point to be inserted.
!
!    Input, integer WP, the index of vertex in PVL which is to be the predecessor
!    vertex of the inserted vertex.
!
!    Input/output, integer NVC, the number of positions used in VCL array.
!
!    Input/output, integer NVERT, the number of positions used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, integer MAXPV, the maximum size available for PVL array.
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Output, integer W, the index of inserted vertex in PVL.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3 or 5.
!
  integer maxpv
  integer maxvc
!
  integer, parameter :: edgv = 4
  double precision iang(maxpv)
  integer ierror
  integer, parameter :: loc = 1
  integer nvc
  integer nvert
  double precision d_pi
  integer, parameter :: polg = 2
  integer pvl(4,maxpv)
  integer, parameter :: succ = 3
  double precision tol
  double precision vcl(2,maxvc)
  integer w
  integer wp
  integer ws
  integer ww
  integer wwp
  integer wws
  double precision xi
  double precision yi
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( nvc >= maxvc ) then
    ierror = 3
    return
  end if

  if ( nvert+2 > maxpv ) then
    ierror = 5
    return
  end if
!
!  Update linked list of vertices of polygon containing vertex WP.
!
  nvc = nvc + 1
  vcl(1,nvc) = xi
  vcl(2,nvc) = yi
  ws = pvl(succ,wp)
  nvert = nvert + 1
  w = nvert
  pvl(loc,w) = nvc
  pvl(polg,w) = pvl(polg,wp)
  pvl(succ,wp) = w
  pvl(succ,w) = ws
  iang(w) = d_pi ( )
  pvl(edgv,w) = pvl(edgv,wp)
!
!  If edge containing (XI,YI) is shared by another polygon,
!  then also update linked list of vertices of that polygon.
!
  if ( pvl(edgv,wp) > 0 ) then
    wws = pvl(edgv,wp)
    wwp = pvl(succ,wws)
    nvert = nvert + 1
    ww = nvert
    pvl(loc,ww) = nvc
    pvl(polg,ww) = pvl(polg,wws)
    pvl(succ,wws) = ww
    pvl(succ,ww) = wwp
    iang(ww) = d_pi ( )
    pvl(edgv,wp) = ww
    pvl(edgv,ww) = wp
    pvl(edgv,wws) = w
  end if

  return
end subroutine




subroutine intpg ( nvrt, xc, yc, ctrx, ctry, arpoly, hflag, umdf, wsq, nev, &
  ifv, listev, ivrt, edgval, vrtval, vcl, mdfint, mean, stdv, mdftr )
!
!*******************************************************************************
!
!! INTPG integrates the mesh distribution function in a convex polygon.
!
!
!  Purpose: 
!
!    Compute integral of MDF2(X,Y) [heuristic mdf] or
!    UMDF(X,Y) [user-supplied mdf] in convex polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices in polygon.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the coordinates of 
!    polygon vertices in counter clockwise order, translated so that 
!    centroid is at origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, double precision CTRX, CTRY, the coordinates of centroid 
!    before translation.
!
!    Input, double precision ARPOLY, the area of polygon.
!
!    Input, logical HFLAG, is .TRUE. if heuristic mdf, .FALSE. if 
!    user-supplied mdf.
!
!    Input, external UMDF, the name of the user supplied MDF routine, of 
!    the form:
!
!      function umdf ( x, y )
!      double precision umdf
!      double precision x
!      double precision y
!
!    Input, double precision WSQ, the square of width of original polygon 
!    of decomposition.
!
!    Input, integer NEV, integer IFV, integer LISTEV(1:NEV), output from 
!    routine PRMDF2.
!
!    Input, IVRT(1:*), double precision EDGVAL(1:*), 
!    double precision VRTVAL(1:*), arrays output from DSMDF2;
!    if .NOT. HFLAG then only first array exists.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, double precision MDFINT, the integral of MDF in polygon.
!
!    Output, double precision MEAN, the mean MDF value in polygon.
!
!    Output, double precision STDV, the standard deviation of MDF in polygon.
!
!    Output, double precision MDFTR(0:NVRT-1), the mean MDF value in each 
!    triangle of polygon;  triangles are determined by polygon vertices 
!    and centroid.
!
  integer nev
  integer, parameter :: nqpt = 3
  integer nvrt
!
  double precision areatr
  double precision arpoly
  double precision ctrx
  double precision ctry
  double precision d
  double precision edgval(*)
  logical hflag
  integer i
  integer ifv
  integer ivrt(*)
  integer j
  integer k
  integer kp1
  integer l
  integer listev(nev)
  integer m
  double precision mdfint
  double precision mdfsqi
  double precision mdftr(0:nvrt-1)
  double precision mean
  double precision, save, dimension ( 3, nqpt ) :: qc = reshape ( (/ &
    0.6666666666666666D+00, 0.1666666666666667D+00, 0.1666666666666667D+00, &
    0.1666666666666667D+00, 0.6666666666666666D+00, 0.1666666666666667D+00, &
    0.1666666666666667D+00, 0.1666666666666667D+00, 0.6666666666666666D+00/), &
    (/ 3, nqpt /) )
  double precision s
  double precision stdv
  double precision sum1
  double precision sum2
  double precision temp
  double precision umdf
  double precision val
  double precision vcl(2,*)
  double precision vrtval(*)
  double precision wsq
  double precision, save, dimension ( nqpt ) :: wt = &
    (/ 0.3333333333333333D+00, 0.3333333333333333D+00, 0.3333333333333333D+00 /)
  double precision x
  double precision x0
  double precision x1
  double precision xc(0:nvrt)
  double precision xx
  double precision y
  double precision y0
  double precision y1
  double precision yc(0:nvrt)
  double precision yy
!
  external umdf
!
!  NQPT is number of quadrature points for numerical integration in triangle.
!  WT(I) is weight of Ith quadrature point.
!  QC(1:3,I) are barycentric coordinates of Ith quadrature point.
!
  mdfint = 0.0D+00
  mdfsqi = 0.0D+00

  do l = 0, nvrt-1

    areatr = 0.5D+00 * ( xc(l) * yc(l+1) - xc(l+1) * yc(l) )
    sum1 = 0.0D+00
    sum2 = 0.0D+00

    do m = 1, nqpt

      xx = qc(1,m) * xc(l) + qc(2,m) * xc(l+1)
      yy = qc(1,m) * yc(l) + qc(2,m) * yc(l+1)
      if ( hflag ) then
!
!             VAL = MDF2(XX+CTRX,YY+CTRY,WSQ,NEV,IFV,LISTEV,IVRT,
!    $            EDGVAL,VRTVAL,VCL)
!
!  Insert code for function MDF2 to reduce number of calls.
!
         x = xx + ctrx
         y = yy + ctry
           s = wsq
           do i = 1, nev
            k = listev(i)
            if ( k < 0 ) then
               k = -k
               d = ( vcl(1,k) - x )**2 + ( vcl(2,k) - y )**2
               d = max ( 0.25D+00 * d, vrtval(k) )
               s = min ( s, d )
            else
               kp1 = k + 1
               if ( i == nev .and. ifv > 0) then
                 kp1 = ifv
               end if
               j = ivrt(kp1)
               x0 = x - vcl(1,j)
               y0 = y - vcl(2,j)
               x1 = vcl(1,ivrt(k)) - vcl(1,j)
               y1 = vcl(2,ivrt(k)) - vcl(2,j)

               if ( x0 * x1 + y0 * y1 <= 0.0D+00 ) then
                 d = x0**2 + y0**2
               else
                 x0 = x0 - x1
                 y0 = y0 - y1
                 if ( x0 * x1 + y0 * y1 >= 0.0D+00 ) then
                   d = x0**2 + y0**2
                 else
                   d = ( x1 * y0 - y1 * x0 )**2 / ( x1**2 + y1**2 )
                 end if
               end if

               d = max ( 0.25D+00 * d, edgval(k) )
               s = min ( s, d )
            end if
           end do
           val = 1.0D+00 / s
      else
        val = umdf ( xx+ctrx, yy+ctry )
      end if

      temp = wt(m) * val
      sum1 = sum1 + temp
      sum2 = sum2 + temp * val

    end do

    mdftr(l) = sum1
    mdfint = mdfint + sum1 * areatr
    mdfsqi = mdfsqi + sum2 * areatr

  end do

  mean = mdfint / arpoly
  stdv = mdfsqi / arpoly - mean**2
  stdv = sqrt ( max ( stdv, 0.0D+00 ) )

  return
end subroutine



subroutine inttri ( nvrt, xc, yc, h, ibot, costh, sinth, ldv, nvc, ntri,  &
  maxvc, maxti, maxcw, vcl, til, ncw, cwalk, ierror )
!
!*******************************************************************************
!
!! INTTRI generates triangles inside a convex polygon.
!
!
!  Purpose: 
!
!    Generate triangles inside convex polygon using quasi-uniform grid of
!    spacing H.  It is assumed that the diameter of the polygon is parallel 
!    to the Y axis.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    in counter clockwise order; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, double precision H, the spacing of mesh vertices in polygon.
!
!    Input, integer IBOT, the index of bottom vertex; diameter contains vertices
!    (XC(0),YC(0)) and (XC(IBOT),YC(IBOT)).
!
!    Input, double precision COSTH, SINTH; COS(THETA), SIN(THETA) where 
!    THETA in [-PI,PI] is rotation angle to get diameter parallel to y-axis.
!
!    Input, integer LDV, the leading dimension of VCL in calling routine.
!
!    Input/output, integer NVC, the number of coordinates or positions 
!    used in VCL array.
!
!    Input/output, integer NTRI, the number of triangles or positions 
!    used in TIL.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, integer MAXTI, the maximum size available for TIL array.
!
!    Input, integer MAXCW, the maximum size available for CWALK array; 
!    assumed to be >= 6*(1 + INT((YC(0) - YC(IBOT))/H)).
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer TIL(1:3,1:NTRI), the triangle incidence list.
!
!    Output, integer NCW, the number of mesh vertices in closed walk, 
!    except NCW = 0 for 1 vertex.
!
!    Output, integer CWALK(0:NCW), indices in VCL of mesh vertices of closed
!    walk; CWALK(0) = CWALK(NCW)
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 9, or 10
!
  integer ldv
  integer maxcw
  integer maxti
  integer maxvc
  integer nvrt
!
  double precision a
  double precision b
  double precision costh
  integer cwalk(0:maxcw)
  double precision cy
  double precision h
  integer i
  integer ibot
  integer ierror
  integer il
  integer im1l
  integer im1r
  integer ir
  integer j
  integer k
  integer l
  integer l0
  integer l1
  integer lw
  integer m
  integer n
  integer ncw
  integer ntri
  integer nvc
  integer p
  integer r
  integer r0
  integer r1
  integer rw
  double precision sinth
  double precision sy
  integer til(3,maxti)
  double precision tol
  double precision vcl(ldv,maxvc)
  double precision x
  double precision xc(0:nvrt)
  double precision xj
  double precision xk
  double precision xl
  double precision xm1l
  double precision xm1r
  double precision xr
  double precision y
  double precision yc(0:nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  n = int ( ( yc(0) - yc(ibot) ) / h )
  y = yc(0) - 0.5D+00 * ( yc(0) - yc(ibot ) - dble ( n ) * h )
  l = 0
  r = nvrt

  do i = 0, n
!
!  Determine left and right x-coordinates of polygon for
!  scan line with y-coordinate Y, and generate mesh vertices.
!
    do while ( yc(l+1) > y )
      l = l + 1
    end do

    do while ( yc(r-1) > y )
      r = r - 1
    end do

    xl = xc(l) + ( xc(l+1) - xc(l) ) * ( y - yc(l) ) / ( yc(l+1) - yc(l) )
    xr = xc(r) + ( xc(r-1) - xc(r) ) * ( y - yc(r) ) / ( yc(r-1) - yc(r) )
    m = int ( ( xr - xl ) / h )
    x = xl + 0.5D+00 * ( xr - xl - dble ( m ) * h )

    if ( nvc + m + 1 > maxvc ) then
      ierror = 3
      return
    end if

    cy = costh * y
    sy = sinth * y
    il = nvc + 1
    xl = x

    do j = 0, m
      nvc = nvc + 1
      vcl(1,nvc) = costh * x + sy
      vcl(2,nvc) = cy - sinth * x
      x = x + h
    end do

    ir = nvc
    xr = x - h

    if ( n == 0 ) then

      ncw = 0
      cwalk(0) = nvc
      return

    else if ( i == 0 ) then

      lw = 0
      cwalk(lw) = il
      rw = maxcw + 1

      do j = il, ir
        rw = rw - 1
        cwalk(rw) = j
      end do

      go to 100

    end if
!
!  Generate triangles between scan lines Y+H and Y.
!
    a = max ( xl, xm1l )
    b = min ( xr, xm1r )

    if ( xm1l == a ) then
      l0 = im1l
      x = ( xm1l - xl ) / h
      j = int(x + tol)
      if ( abs ( x - dble ( j ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      l1 = il + j
    else
      l1 = il
      x = ( xl - xm1l ) / h
      j = int ( x + tol )
      if ( abs ( x - dble ( j ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      l0 = im1l + j
    end if

    if ( xm1r == b ) then
      r0 = im1r
      x = ( xr - xm1r ) / h
      j = int ( x + tol )
      if ( abs ( x - dble ( j ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      r1 = ir - j
    else
      r1 = ir
      x = ( xm1r - xr ) / h
      j = int ( x + tol )
      if ( abs ( x - dble(j) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      r0 = im1r - j
    end if

    if ( l0 < r0 .or. l1 < r1 ) then

      j = l0
      k = l1
      xj = xm1l + dble ( j-im1l ) * h
      xk = xl + dble ( k - il ) * h

      do

        if ( k < r1 .and. ( xk <= xj .or. j == r0 ) ) then
          p = k
          k = k + 1
          xk = xk + h
        else
          p = j
          j = j + 1
          xj = xj + h
        end if

        ntri = ntri + 1

        if ( ntri > maxti ) then
          ierror = 9
          return
        end if

        til(1,ntri) = j
        til(2,ntri) = p
        til(3,ntri) = k

        if ( j >= r0 .and. k >= r1 ) then
          exit
        end if

      end do

    end if
!
!  Generate paths of closed walk between scan lines Y+H and Y.
!
    if ( xm1l < xl ) then
      do j = im1l+1, l0
        lw = lw + 1
        cwalk(lw) = j
      end do
      lw = lw + 1
      cwalk(lw) = il
    else
      do j = l1, il, -1
        lw = lw + 1
        cwalk(lw) = j
      end do
    end if

    if ( xm1r > xr ) then
      do j = im1r-1, r0, -1
        rw = rw - 1
        cwalk(rw) = j
      end do
      rw = rw - 1
      cwalk(rw) = ir
    else
      do j = r1, ir
        rw = rw - 1
        cwalk(rw) = j
      end do
    end if

100 continue

    y = y - h
    im1l = il
    im1r = ir
    xm1l = xl
    xm1r = xr

  end do
!
!  Add last path of left walk and shift indices of right walk.
!
  if ( m == 0 ) then
    rw = rw + 1
  else
    do j = il+1, ir-1
      lw = lw + 1
      cwalk(lw) = j
    end do
  end if

  if ( rw <= lw ) then
    ierror = 10
    return
  end if

  do j = rw, maxcw
    lw = lw + 1
    cwalk(lw) = cwalk(j)
  end do

  ncw = lw

  return
end subroutine



subroutine isftdw ( l, u, k, lda, a, map )
!
!*******************************************************************************
!
!! ISFTDW sifts A(*,MAP(L)) down a heap of size U.
!
!
!  Purpose: 
!
!    Sift A(*,MAP(L)) down a heap of size U.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer L, U, the lower and upper index of part of heap.
!
!    Input, integer K, the dimension of points.
!
!    Input, integer LDA, the leading dimension of array A in calling routine.
!
!    Input, integer A(1:K,1:*); see routine IHPSRT.
!
!    Input/output, integer MAP(1:*); see routine IHPSRT.
!
  integer lda
!
  integer a(lda,*)
  integer i
  logical iless
  integer j
  integer k
  integer l
  integer map(*)
  integer u
  integer t
!
  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( iless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( iless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end subroutine



subroutine ivec_identity ( n, a )
!
!*******************************************************************************
!
!! IVEC_IDENTITY sets an integer vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 November 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n
    a(i) = i
  end do

  return
end subroutine



subroutine jnhole ( itophv, angspc, angtol, nvc, nvert, maxvc, maxpv,  &
  maxiw, maxwk, vcl, hvl, pvl, iang, iwk, wk, ierror )
!
!*******************************************************************************
!
!! JNHOLE joins a hole boundary to the boundary of a polygon.
!
!
!  Purpose: 
!
!    Join hole boundary to boundary of polygon containing hole
!    by finding a cut edge originating from the top vertex of hole
!    to a point on outer polygon boundary above it.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer ITOPHV, the index in PVL of top vertex of hole.
!
!    Input, double precision ANGSPC, the angle spacing parameter used in 
!    controlling the vertices to be considered as an endpoint of a separator.
!
!    Input, double precision ANGTOL, the angle tolerance parameter used 
!    in accepting separator(s).
!
!    Input/output, integer NVC, the number of positions used in VCL array.
!
!    Input/output, integer NVERT, the number of positions used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, integer MAXPV, the maximum size available for PVL array.
!
!    Input, integer MAXIW, the maximum size available for IWK array; should 
!    be about 3 times number of vertices in outer polygon.
!
!    Input, integer MAXWK, the maximum size available for WK array; should 
!    be about 5 times number of vertices in outer polygon.
!
!    Input/output, double precision VCL(1:2,1:NVC),the vertex coordinate list.
!
!    Input, integer HVL(1:*), the head vertex list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT), 
!    the polygon vertex list and interior angles.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 5, 6, 7, 206 to 210, 212, 218, or 219.
!
  integer maxiw
  integer maxpv
  integer maxvc
  integer maxwk
!
  double precision angle
  double precision angspc
  double precision angtol
  double precision dy
  integer, parameter :: edgv = 4
  integer hv
  integer hvl(*)
  double precision iang(maxpv)
  integer ierror
  integer ilft
  integer ipoly
  integer irgt
  integer itophv
  integer iv
  integer ivs
  integer iwk(maxiw)
  integer l
  integer, parameter :: loc = 1
  integer lv
  integer lw
  integer, parameter :: msglvl = 0
  integer nvc
  integer nvert
  double precision d_pi
  integer, parameter :: polg = 2
  integer pvl(4,maxpv)
  double precision s
  double precision slft
  double precision srgt
  integer, parameter :: succ = 3
  integer succil
  integer succir
  double precision tol
  double precision vcl(2,maxvc)
  integer vp
  integer vr
  integer vs
  integer vv
  integer w
  double precision wk(maxwk)
  integer ww
  double precision xint
  double precision xlft
  double precision xrgt
  double precision xt
  double precision xv
  double precision xvs
  double precision ylft
  double precision yrgt
  double precision yt
  double precision ytmtol
  double precision ytptol
  double precision yv
  double precision yvs
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( nvc+3 > maxvc ) then
    ierror = 3
    return
  end if

  if ( nvert+5 > maxpv ) then
    ierror = 5
    return
  end if
!
!  Determine 'closest' vertices on outer boundary which are to the
!  left and right of top vertex of hole and on the horizontal line
!  through top vertex. The two closest vertices must be on edges
!  which intersect the horizontal line and are partially above the
!  line. Ties are broken (in the case of a vertex on a cut edge)
!  by choosing the vertex on the edge of maximum or minimum dx/dy
!  slope depending on whether the vertex is to the left or right
!  of top vertex, respectively.
!
  ipoly = pvl(polg,itophv)
  lv = pvl(loc,itophv)
  xt = vcl(1,lv)
  yt = vcl(2,lv)
  dy = 0.0D+00
  hv = hvl(ipoly)
  iv = hv
  yv = vcl(2,pvl(loc,iv))

  do

    iv = pvl(succ,iv)
    yvs = vcl(2,pvl(loc,iv))
    dy = max ( dy, abs ( yvs - yv ) )
    yv = yvs

    if ( iv == hv ) then
      exit
    end if

  end do

  ytmtol = yt - tol * dy
  ytptol = yt + tol * dy
  ilft = 0
  irgt = 0
  xlft = 0.0D+00
  xrgt = 0.0D+00
  hv = hvl(ipoly)
  iv = hv
  lv = pvl(loc,iv)
  xv = vcl(1,lv)
  yv = vcl(2,lv)

20 continue

  ivs = pvl(succ,iv)
  lv = pvl(loc,ivs)
  xvs = vcl(1,lv)
  yvs = vcl(2,lv)

  if ( yv <= ytptol .and. yvs > ytptol ) then

    if ( yv >= ytmtol ) then
      if ( xv > xt ) then
        if ( xv < xrgt .or. irgt == 0 ) then
          irgt = iv
          xrgt = xv
          yrgt = yv
          srgt = (xvs - xv ) / ( yvs - yv )
        else if ( xv == xrgt ) then
          s = ( xvs - xv ) / ( yvs - yv )
          if ( s < srgt ) then
            irgt = iv
            yrgt = yv
            srgt = s
          end if
        end if
      end if
    else
      xint = ( yt - yv ) * ( xvs - xv ) / ( yvs - yv ) + xv
      if ( xint > xt ) then
        if ( xint < xrgt .or. irgt == 0 ) then
          irgt = iv
          xrgt = xint
          yrgt = yt
        end if
      end if
    end if

  else if ( yv > ytptol .and. yvs <= ytptol ) then

    if ( yvs >= ytmtol ) then
      if ( xvs < xt ) then
        if ( xvs > xlft .or. ilft == 0 ) then
          ilft = iv
          xlft = xvs
          ylft = yvs
          slft = ( xvs - xv ) / ( yvs - yv )
        else if ( xvs == xlft ) then
          s = ( xvs - xv ) / ( yvs - yv )
          if ( s > slft ) then
            ilft = iv
            ylft = yvs
            slft = s
          end if
        end if
      end if
    else
      xint = ( yt - yv ) * ( xvs - xv ) / ( yvs - yv ) + xv
      if ( xint < xt ) then
        if ( xint > xlft .or. ilft == 0 ) then
          ilft = iv
          xlft = xint
          ylft = yt
        end if
      end if
    end if
  end if

  iv = ivs
  xv = xvs
  yv = yvs

  if ( iv /= hv ) then
    go to 20
  end if

  if ( ilft == 0 .or. irgt  ==  0 ) then
   ierror = 218
   return
  end if
!
!  Temporarily modify PVL to pass the subregion 'above' top vertex
!  of hole to routine RESVRT. The top vertex is the reflex vertex
!  passed to RESVRT (in the temporary subregion, it has interior
!  angle PI). This causes one separator to be chosen by RESVRT
!  and its other endpoint is above the top vertex.
!
  succil = pvl(succ,ilft)
  succir = pvl(succ,irgt)
  vcl(1,nvc+2) = xlft
  vcl(2,nvc+2) = ylft
  vcl(1,nvc+3) = xrgt
  vcl(2,nvc+3) = yrgt
  vp = nvert + 3
  vr = nvert + 4
  vs = nvert + 5

  iang(vr) = angle ( xlft, ylft, xt, yt, xrgt, yrgt )

  if ( iang(vr) < d_pi ( ) - tol .or. iang(vr) > d_pi ( ) + tol ) then
    ierror = 219
    return
  end if

  pvl(loc,vp) = nvc + 2
  pvl(polg,vp) = ipoly
  pvl(succ,vp) = vr
  pvl(edgv,vp) = 0
  pvl(loc,vr) = pvl(loc,itophv)
  pvl(polg,vr) = ipoly
  pvl(succ,vr) = vs
  pvl(edgv,vr) = 0
  pvl(loc,vs) = nvc + 3
  pvl(polg,vs) = ipoly
  pvl(succ,vs) = succir
  pvl(edgv,vs) = pvl(edgv,irgt)
  pvl(succ,ilft) = vp
  lv = pvl(loc,ilft)
  iang(vp) = angle ( vcl(1,lv), vcl(2,lv), xlft, ylft, xt, yt )
  lv = pvl(loc,succir)
  iang(vs) = angle ( xt, yt, xrgt, yrgt, vcl(1,lv), vcl(2,lv) )
  w = 0

  call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, maxwk, &
    vcl, pvl, iang, w, ww, iwk, wk, ierror )
!
!  Remove temporary modification to PVL. There are three cases
!  depending on where the endpoint of separator is located:
!  successor of closest vertex to the right of top vertex,
!  predecessor of closest vertex to the left of top vertex,
!  or neither of these.
!
  if ( pvl(succ,vs) == w ) then
    pvl(succ,ilft) = succil
    pvl(succ,irgt) = w
    pvl(edgv,irgt) = pvl(edgv,vs)
    if ( pvl(edgv,irgt) > 0 ) then
      pvl(edgv,pvl(edgv,irgt)) = irgt
    end if
  else if ( pvl(succ,ilft) == w ) then
    pvl(succ,w) = succil
  else
    pvl(succ,ilft) = succil
  end if

  if ( ierror /= 0 ) then
    return
  end if
!
!  Update PVL with cut edge, i.e. join linked lists of vertices
!  of the hole polygon and the outer boundary polygon into one
!  linked list of vertices by adding the cut edge from the top
!  vertex of hole to the vertex on the outer boundary.
!
  nvert = nvert + 2
  vv = nvert - 1
  ww = nvert
  lv = pvl(loc,itophv)
  lw = pvl(loc,w)
  pvl(loc,vv) = lv
  pvl(loc,ww) = lw
  pvl(polg,vv) = ipoly
  pvl(polg,ww) = ipoly
  pvl(succ,vv) = pvl(succ,itophv)
  pvl(succ,ww) = pvl(succ,w)
  pvl(succ,itophv) = ww
  pvl(succ,w) = vv
  pvl(edgv,vv) = pvl(edgv,itophv)
  pvl(edgv,ww) = pvl(edgv,w)
  pvl(edgv,itophv) = w
  pvl(edgv,w) = itophv

  if ( pvl(edgv,vv) > 0 ) then
    pvl(edgv,pvl(edgv,vv)) = vv
  end if

  if ( pvl(edgv,ww) > 0 ) then
    pvl(edgv,pvl(edgv,ww)) = ww
  end if

  l = pvl(loc,pvl(succ,vv))
  iang(vv) = angle ( vcl(1,lw), vcl(2,lw), vcl(1,lv), vcl(2,lv), vcl(1,l), &
    vcl(2,l) )
  iang(itophv) = iang(itophv) - iang(vv)
  l = pvl(loc,pvl(succ,ww))
  iang(ww) = angle ( vcl(1,lv), vcl(2,lv), vcl(1,lw), vcl(2,lw), vcl(1,l), &
    vcl(2,l) )
  iang(w) = iang(w) - iang(ww)

  if ( msglvl == 2 ) then
    write ( *,600) itophv,w,vcl(1,lv),vcl(2,lv),vcl(1,lw),vcl(2,lw)
  end if

600 format (1x,2i7,4f15.7)

  return
end subroutine


subroutine lop ( itr, ind, mxedg, top, ldv, vcl, til, tedg, sptr )
!
!*******************************************************************************
!
!! LOP applies the local optimization procedure to two triangles.
!
!
!  Purpose: 
!
!    Apply local optimization procedure to two triangles
!    indicated by ITR(1) and ITR(2). This may result in swapping
!    diagonal edge of quadrilateral.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer ITR(1:2), the indices of triangles for LOP.
!
!    Input, integer IND(1:2), indices indicating common edge of triangles.
!
!    Input, integer MXEDG, the maximum index of edge to be considered for LOP.
!
!    Input/output, integer TOP, the index of SPTR indicating top of stack.
!
!    Input, integer LDV, the leading dimension of VCL in calling routine.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer TIL(1:3,1:*), the triangle incidence list.
!
!    Input/output, integer TEDG(1:3,1:*), the triangle edge indices; 
!    see routine CVDTRI.
!
!    Input/output, integer SPTR(1:*), stack pointers; see routine CVDTRI.
!
  integer ldv
!
  integer a
  integer b
  integer c
  integer d
  integer diaedg
  integer i
  integer i_wrap
  integer iedg
  integer in
  integer ind(2)
  integer ind1m1
  integer ind1p1
  integer ind2m1
  integer ind2p1
  integer itr(2)
  integer j
  integer mxedg
  integer top
  integer sptr(*)
  integer tedg(3,*)
  integer til(3,*)
  double precision vcl(ldv,*)
!
!  Common edge is BC, other two vertices are A and D.
!
  iedg = tedg(ind(1),itr(1))
  sptr(iedg) = -1

  ind1m1 = i_wrap ( ind(1) - 1, 1, 3 )
  ind1p1 = i_wrap ( ind(1) + 1, 1, 3 )
  ind2m1 = i_wrap ( ind(2) - 1, 1, 3 )
  ind2p1 = i_wrap ( ind(2) + 1, 1, 3 )

  b = til(ind(1),itr(1))
  c = til(ind1p1,itr(1))
  a = til(ind1m1,itr(1))
  d = til(ind2m1,itr(2))

  in = diaedg ( vcl(1,d), vcl(2,d), vcl(1,c), vcl(2,c), vcl(1,a), vcl(2,a), &
    vcl(1,b), vcl(2,b) )

  if ( in == 1 ) then
!
!  Check if four edges of quadrilateral should be put on LOP
!  stack, and swap edge BC for AD.
!
   i = tedg(ind1m1,itr(1))

   do j = 1, 4

      if ( j == 2 ) then
        i = tedg(ind1p1,itr(1))
      else if ( j == 3 ) then
        i = tedg(ind2m1,itr(2))
      else if ( j == 4 ) then
        i = tedg(ind2p1,itr(2))
      end if

      if ( i <= mxedg ) then
        if ( sptr(i) == -1 ) then
          sptr(i) = top
          top = i
        end if
      end if

    end do

    til(ind1p1,itr(1)) = d
    til(ind2p1,itr(2)) = a
    tedg(ind(1),itr(1)) = tedg(ind2p1,itr(2))
    tedg(ind(2),itr(2)) = tedg(ind1p1,itr(1))
    tedg(ind1p1,itr(1)) = iedg
    tedg(ind2p1,itr(2)) = iedg

  end if

  return
end subroutine


function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )
!
!*******************************************************************************
!
! LRLINE determines if a point is left of, right or, or on a directed line.
!
!
!  Discussion:
!
!    The directed line is paralled to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, double precision XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, double precision DV, the signed distance of the directed line 
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).  
!    DV is positive for a line to the left of the base line.
! 
!    Output, integer LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  double precision dv
  double precision dx
  double precision dxu
  double precision dy
  double precision dyu
  integer lrline
  double precision t
  double precision tol
  double precision tolabs
  double precision xu
  double precision xv1
  double precision xv2
  double precision yu
  double precision yv1
  double precision yv2
!
  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end function

subroutine lufac ( a, lda, n, tol, ipvt, singlr )
!
!*******************************************************************************
!
!! LUFAC computes the LU factorization of a matrix.
!
!
!  Purpose: 
!
!    Obtain LU factorization of matrix A, i.e. apply Gaussian
!    elimination with partial pivoting to A.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input/output, double precision A(1:N,1:N), the matrix.  On input, 
!    the N by N matrix to be factored.  On output, the upper triangular 
!    matrix U and multipliers of unit lower triangular matrix L (if matrix 
!    A is nonsingular).
!
!    Input, integer LDA, the leading dimension of array A in calling routine.
!
!    Input, integer N, the order of matrix A.
!
!    Input, double preicison TOL, the relative tolerance for detecting 
!    singularity of A.
!
!    Output, integer IPVT(1:N-1), the pivot indices.
!
!    Output, logical SINGLR, is .TRUE. if matrix is singular; this occurs 
!    when the magnitude of a pivot element is <= TOL*MAX(|A(I,J)|).
!
  integer lda
  integer n
!
  double precision a(lda,n)
  integer i
  integer ipvt(n-1)
  integer j
  integer k
  integer kp1
  integer m
  logical singlr
  double precision t
  double precision tol
  double precision tolabs
!
  if ( n < 1 ) then
    return
  end if

  singlr = .true.

  t = 0.0D+00
  do j = 1, n
    do i = 1, n
      t = max ( t, abs ( a(i,j) ) )
    end do
  end do

  tolabs = tol * t

  do k = 1, n-1

    kp1 = k + 1
    m = k

    do i = k+1, n
      if ( abs ( a(i,k)) > abs ( a(m,k) ) ) then
        m = i
      end if
    end do

    ipvt(k) = m
    t = a(m,k)
    a(m,k) = a(k,k)
    a(k,k) = t

    if ( abs ( t) <= tolabs ) then
      return
    end if

    do i = k+1, n
      a(i,k) = a(i,k) / t
    end do

    do j = k+1, n

      t = a(m,j)
      a(m,j) = a(k,j)
      a(k,j) = t

      if ( t /= 0.0E+00 ) then

        do i = k+1, n
          a(i,j) = a(i,j) - a(i,k) * t
        end do

      end if

    end do

  end do

  if ( abs ( a(n,n) ) > tolabs ) then
    singlr = .false.
  end if

  return
end
subroutine lusol ( a, lda, n, ipvt, b )
!
!*******************************************************************************
!
!! LUSOL solves a linear system with an LU factored matrix.
!
!
!  Purpose: 
!
!    Solve linear system A*X = B given LU factorization of A.
!    It is assumed that A is nonsingular.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision A(1:N,1:N), contains factors L, U output 
!    from routine LUFAC.
!
!    Input, integer LDA, the leading dimension of array A in calling routine.
!
!    Input, integer N, the order of matrix A.
!
!    Input, integer IPVT(1:N-1), the pivot indices from routine LUFAC.
!
!    Input/output, B(1:N).  On input, the right hand side vector.
!    On output, the solution vector X.
!
  integer lda
  integer n
!
  double precision a(lda,n)
  double precision b(n)
  integer ipvt(n-1)
  integer i
  integer k
  integer m
  double precision t
!
!  Forward elimination
!
  do k = 1, n-1

    m = ipvt(k)
    t = b(m)
    b(m) = b(k)
    b(k) = t

    do i = k+1, n
      b(i) = b(i) - a(i,k) * t
    end do

  end do
!
!  Back substitution
!
  do k = n, 2, -1

    t = b(k) / a(k,k)
    b(k) = t

    do i = 1, k-1
      b(i) = b(i) - a(i,k) * t
    end do

  end do

  b(1) = b(1) / a(1,1)

  return
end
function mdf2 ( x, y, wsq, nev, ifv, listev, ivrt, edgval, vrtval, vcl )
!
!*******************************************************************************
!
!! MDF2 evaluates the heuristic mesh distribution function at (X,Y).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of point.
!
!    Input, double precision WSQ, the square of the width of the polygon 
!    containing (X,Y).
!
!    Input, integer NEV, IFV, LISTEV(1:NEV), output from routine PRMDF2.
!
!    Input, integer IVRT(1:*), output from DSMDF2.
!
!    Input, double precision EDGVAL(1:*), output from DSMDF2.
!
!    Input, double precision VRTVAL(1:*), output from DSMDF2.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, double precision MDF2, the reciprocal of square of length 
!    scale at (X,Y).
!
  integer nev
!
  double precision d
  double precision edgval(*)
  integer i
  integer ifv
  integer ivrt(*)
  integer j
  integer k
  integer kp1
  integer listev(nev)
  double precision mdf2
  double precision s
  double precision vcl(2,*)
  double precision vrtval(*)
  double precision wsq
  double precision x
  double precision x0
  double precision x1
  double precision y
  double precision y0
  double precision y1
!
  s = wsq

  do i = 1, nev

    k = listev(i)

    if ( k < 0 ) then
      k = -k
      d = ( vcl(1,k) - x )**2 + ( vcl(2,k) - y )**2
      d = max ( 0.25D+00 * d, vrtval(k) )
      s = min ( s, d )
    else
      kp1 = k + 1
      if ( i == nev .and. ifv > 0 ) then
        kp1 = ifv
      end if
      j = ivrt(kp1)
      x0 = x - vcl(1,j)
      y0 = y - vcl(2,j)
      x1 = vcl(1,ivrt(k)) - vcl(1,j)
      y1 = vcl(2,ivrt(k)) - vcl(2,j)

      if ( x0 * x1 + y0 * y1 <= 0.0D+00 ) then
        d = x0**2 + y0**2
      else
        x0 = x0 - x1
        y0 = y0 - y1
        if ( x0 * x1 + y0 * y1 >= 0.0D+00 ) then
          d = x0**2 + y0**2
        else
          d = ( x1 * y0 - y1 * x0 )**2 / ( x1**2 + y1**2 )
        end if
      end if

      d = max ( 0.25D+00 * d, edgval(k) )
      s = min ( s, d )

    end if

  end do

  mdf2 = 1.0D+00 / s

  return
end
subroutine mfdec2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid,  &
  nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl,  &
  pvl, iang, ivrt, xivrt, widsq, edgval, vrtval, area, psi, iwk, wk, ierror )
!
!*******************************************************************************
!
!! MFDEC2 subdivides polygons to decrease mesh distribution variation.
!
!
!  Purpose: 
!
!    Further subdivide convex polygons so that the variation
!    of heuristic or user-supplied mesh distribution function in
!    each polygon is limited.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical HFLAG, is .TRUE. if heuristic mdf, .FALSE. if 
!    user-supplied mdf.
!
!    Input, external UMDF, a user-supplied mdf, of the form
!
!      function umdf ( x, y )
!      double precision umdf
!      double precision x
!      double precision y
!
!    Input, double precision KAPPA, the mesh smoothness parameter in 
!    interval [0.0,1.0].
!
!    Input, double precision ANGSPC, the angle spacing parameter in radians 
!    used to determine extra points as possible endpoints of separators.
!
!    Input, double precision ANGTOL, the angle tolerance parameter in 
!    radians used in accepting separators.
!
!    Input, double precision DMIN, a parameter used to determine if 
!    variation of mdf in polygon is 'sufficiently high'.
!
!    Input, integer NMIN, a parameter used to determine if 'sufficiently large'
!    number of triangles in polygon.
!
!    Input, integer NTRID, the desired number of triangles in mesh.
!
!    Input/output, integer NVC, the number of vertex coordinates or 
!    positions used in VCL array.
!
!    Input/output, integer NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, integer MAXHV, the maximum size available for HVL, REGNUM, AREA, 
!    PSI arrays.
!
!    Input, integer MAXPV, the maximum size available for PVL, IANG arrays.
!
!    Input, integer MAXIW, the maximum size available for IWK array; should 
!    be about 3*NVRT + INT(2*PI/ANGSPC) where NVRT is maximum number of
!    vertices in a convex polygon of the (input) decomposition.
!
!    Input, integer MAXWK, the maximum size available for WK array; should 
!    be about NPOLG + 3*(NVRT + INT(2*PI/ANGSPC)) + 2.
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, HVL(1:NPOLG), the head vertex list.
!
!    Input/output, PVL(1:4,1:NVERT), IANG(1:NVERT), the polygon vertex list and
!    interior angles.
!
!    Input, integer IVRT(1:NVERT), integer XIVRT(1:NPOLG+1), 
!    double precision WIDSQ(1:NPOLG), double precision EDGVAL(1:NVERT),
!    double precision VRTVAL(1:NVC), arrays output from routine DSMDF2;
!    if .NOT. HFLAG then only first two arrays exist.
!
!    Input/output, double precision AREA(1:NPOLG), the area of convex polygons 
!    in decomposition.
!
!    Output, double precision PSI(1:NPOLG), the mean mdf values in the
!    convex polygons
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 200, or 222.
!
  integer maxhv
  integer maxiw
  integer maxpv
  integer maxvc
  integer maxwk
  integer nvc
  integer nvert
!
  double precision alpha
  double precision angsp2
  double precision angspc
  double precision angtol
  double precision area(maxhv)
  double precision areapg
  double precision arearg
  double precision c1
  double precision c2
  double precision cosalp
  double precision ctrx
  double precision ctry
  double precision delta
  double precision dmin
  double precision dx
  double precision dy
  double precision edgval(nvert)
  logical hflag
  integer hvl(maxhv)
  integer i
  integer i1
  integer i2
  double precision iang(maxpv)
  integer ierror
  integer ifv
  integer ii
  integer inc
  integer indpvl
  double precision intreg
  integer ivrt(nvert)
  integer iwk(maxiw)
  integer j
  integer k
  double precision kappa
  integer l
  integer listev
  integer, parameter :: loc = 1
  integer m
  integer maxn
  double precision mdfint
  integer mdftr
  double precision mean
  integer nev
  integer nmin
  integer np
  integer npolg
  integer ntrid
  double precision numer
  integer nvrt
  double precision nwarea
  integer p
  double precision d_pi
  double precision pi2
  double precision psi(maxhv)
  integer pvl(4,maxpv)
  double precision r
  integer regnum(maxhv)
  double precision sinalp
  double precision stdv
  integer, parameter :: succ = 3
  double precision sumx
  double precision sumy
  double precision theta1
  double precision theta2
  double precision tol
  double precision umdf
  integer v
  double precision vcl(2,maxvc)
  double precision vrtval(nvc)
  integer w
  double precision widsq(npolg)
  double precision wk(maxwk)
  double precision wsq
  double precision x1
  double precision x2
  integer xc
  integer xivrt(npolg+1)
  double precision y1
  double precision y2
  integer yc
!
  external umdf
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  WK(1:NPOLG) is used for MDF standard deviation in polygons.
!  Compute AREARG = area of region and INTREG = estimated integral
!  of MDF2(X,Y) or UMDF(X,Y).
!
  nvrt = 0
  do i = 1, npolg
    nvrt = max ( nvrt, xivrt(i+1) - xivrt(i) )
  end do

  if ( hflag .and. 2 * nvrt > maxiw ) then
    ierror = 6
    return
  end if

  if ( npolg + 3 * nvrt + 2 > maxwk ) then
    ierror = 7
    return
  end if

  listev = 1
  xc = npolg + 1
  yc = xc + nvrt + 1
  mdftr = yc + nvrt + 1
  arearg = 0.0D+00
  intreg = 0.0D+00
  nev = -1

  do i = 1, npolg

    if ( hflag ) then
      wsq = widsq(i)
      call prmdf2 ( i, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, iwk(listev) )
    end if

    if ( nev == 0 ) then
      psi(i) = 1.0D+00 / wsq
      wk(i) = 0.0D+00
      mdfint = psi(i) * area(i)
    else

      nvrt = xivrt(i+1) - xivrt(i)
      k = xivrt(i)

      sumx = 0.0D+00
      sumy = 0.0D+00
      do j = 0, nvrt-1
        l = ivrt(k)
        wk(xc+j) = vcl(1,l)
        wk(yc+j) = vcl(2,l)
        sumx = sumx + wk(xc+j)
        sumy = sumy + wk(yc+j)
        k = k + 1
      end do

      ctrx = sumx / dble ( nvrt )
      ctry = sumy / dble ( nvrt )

      do j = 0, nvrt-1
        wk(xc+j) = wk(xc+j) - ctrx
        wk(yc+j) = wk(yc+j) - ctry
      end do

      wk(xc+nvrt) = wk(xc)
      wk(yc+nvrt) = wk(yc)

      call intpg ( nvrt, wk(xc), wk(yc), ctrx, ctry, area(i), hflag, umdf, &
        wsq, nev, ifv, iwk(listev), ivrt, edgval, vrtval, vcl, mdfint, &
        psi(i), wk(i), wk(mdftr) )

    end if

    arearg = arearg + area(i)
    intreg = intreg + mdfint

  end do
!
!  If HFLAG, compute mean mdf values from KAPPA, etc. Scale PSI(I)'s
!  so that integral in region is 1. Determine which polygons need to
!  be further subdivided (indicated by negative PSI(I) value).
!
  if ( hflag ) then
    c1 = ( 1.0D+00 - kappa ) / intreg
    c2 = kappa / arearg
  else
    c1 = 1.0D+00 / intreg
    c2 = 0.0D+00
  end if

  do i = 1, npolg

    psi(i) = psi(i) * c1 + c2

    if ( c1 * wk(i) > psi(i) * dmin ) then
      if ( ntrid * psi(i) * area(i) > nmin ) then
        psi(i) = -psi(i)
      end if
    end if

  end do
!
!  Further subdivide polygons for which STDV/MEAN > DMIN and
!  (estimated number of triangles) > NMIN.
!
  angsp2 = 2.0D+00 * angspc
  pi2 = 2.0D+00 * d_pi ( )
  inc = int ( pi2 / angspc )
  nev = 0
  np = npolg
  xc = 1

  do i = 1, np

    if ( psi(i) < 0.0D+00 ) then

      if ( hflag ) then
        wsq = widsq(i)
        call prmdf2 ( i, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, &
          iwk(listev) )

      end if

      l = npolg + 1
      k = i

60    continue

      if ( k > npolg ) then
        go to 130
      end if

70       continue

         if ( psi(k) >= 0.0D+00 ) go to 120

          nvrt = 0
          sumx = 0.0D+00
          sumy = 0.0D+00
          j = hvl(k)

80           continue

             nvrt = nvrt + 1
             m = pvl(loc,j)
             sumx = sumx + vcl(1,m)
             sumy = sumy + vcl(2,m)
             j = pvl(succ,j)
          if ( j /= hvl(k) ) go to 80

          ctrx = sumx / dble(nvrt)
          ctry = sumy / dble(nvrt)
          maxn = nvrt + inc

          if ( nev + maxn + 1 > maxiw ) then
            ierror = 6
            return
          else if ( 3*maxn + 2 > maxwk ) then
            ierror = 7
            return
          end if

          yc = xc + maxn + 1
          mdftr = yc + maxn + 1
          indpvl = listev + nev
          nvrt = 0
          m = pvl(loc,j)
          x1 = vcl(1,m) - ctrx
          y1 = vcl(2,m) - ctry
          wk(xc) = x1
          wk(yc) = y1
          theta1 = atan2(y1,x1)
          p = j
          iwk(indpvl) = j

90        continue

             j = pvl(succ,j)
             m = pvl(loc,j)
             x2 = vcl(1,m) - ctrx
             y2 = vcl(2,m) - ctry

             theta2 = atan2 ( y2, x2 )
             if ( theta2 < theta1 ) then
               theta2 = theta2 + pi2
             end if

             delta = theta2 - theta1

             if ( delta >= angsp2 ) then

               m = int ( delta / angspc )
               delta = delta / dble ( m )
               dx = x2 - x1
               dy = y2 - y1
               numer = x1 * dy - y1 * dx
               alpha = theta1

               do ii = 1, m-1
                 alpha = alpha + delta
                 cosalp = cos(alpha)
                 sinalp = sin(alpha)
                 r = numer / ( dy * cosalp - dx * sinalp )
                 nvrt = nvrt + 1
                 wk(xc+nvrt) = r * cosalp
                 wk(yc+nvrt) = r * sinalp
                 iwk(indpvl+nvrt) = -p
               end do

             end if

             nvrt = nvrt + 1
             wk(xc+nvrt) = x2
             wk(yc+nvrt) = y2
             x1 = x2
             y1 = y2
             theta1 = theta2
             p = j
             iwk(indpvl+nvrt) = j

          if ( j /= hvl(k) ) then
            go to 90
          end if


            call intpg ( nvrt, wk(xc), wk(yc), ctrx, ctry, area(k), hflag, &
              umdf, wsq, nev, ifv, iwk(listev), ivrt, edgval, vrtval, &
              vcl, mdfint, mean, stdv, wk(mdftr) )

            psi(k) = mean * c1 + c2

            if ( c1 * stdv > psi(k) * dmin ) then

               if ( ntrid * psi(k) * area(k) > nmin ) then

              call sepmdf ( angtol, nvrt, wk(xc), wk(yc), area(k), &
                          mean, wk(mdftr), iwk(indpvl), iang, i1, i2 )

              if ( i1 < 0 ) then

                if ( yc + 3*nvrt > maxwk ) then
                  ierror = 7
                  return
                end if

                call sepshp ( angtol, nvrt, wk(xc), wk(yc), &
                  iwk(indpvl), iang, i1, i2, wk(yc+nvrt+1), ierror )

                if ( ierror /= 0 ) then
                  return
                end if

              end if

              if ( i1 < 0 ) then
                ierror = 222
                return
              end if

              v = iwk(indpvl+i1)

              if ( v < 0 ) then
                 call insvr2 ( wk(xc+i1)+ctrx, wk(yc+i1)+ctry, -v, &
                   nvc, nvert, maxvc, maxpv, vcl, pvl, iang, v, ierror )
                 if ( ierror /= 0 ) then
                   return
                 end if
              end if

              w = iwk(indpvl+i2)

              if ( w < 0 ) then
                 call insvr2 ( wk(xc+i2)+ctrx, wk(yc+i2)+ctry, -w, &
                   nvc, nvert, maxvc, maxpv, vcl, pvl, iang, w, ierror )
                 if ( ierror /= 0 ) then
                   return
                 end if
              end if

          call insed2 ( v, w, npolg, nvert, maxhv, maxpv, vcl, &
            regnum, hvl, pvl, iang, ierror )

          if ( ierror /= 0 ) then
            return
          end if

          nvrt = 0
          j = hvl(k)

          do

            m = pvl(loc,j)
            wk(xc+nvrt) = vcl(1,m)
            wk(yc+nvrt) = vcl(2,m)
            nvrt = nvrt + 1
            j = pvl(succ,j)

            if ( j == hvl(k) ) then
              exit
            end if

          end do

          nwarea = areapg ( nvrt, wk(xc), wk(yc) ) * 0.5D+00
          area(npolg) = area(k) - nwarea
          area(k) = nwarea
          psi(k) = -psi(k)
          psi(npolg) = psi(k)

        end if

      end if

      go to 70

120   continue

      if ( k == i ) then
        k = l
      else
        k = k + 1
      end if

      go to 60

130   continue

    end if

  end do

  return
end
function minang ( xr, yr, xs, ys, ind, alpha, theta, vcl, pvl, iang )
!
!*******************************************************************************
!
!! MINANG determines the minimum of the boundary angles for a separator.
!
!
!  Purpose: 
!
!    Determine the minimum of the 4 angles at the boundary
!    resulting from using edge joining vertices (XR,YR) and
!    (XS,YS) as a separator.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XR, YR, the coordinates of the reflex vertex.
!
!    Input, double precision XS, YS, the coordinates of other endpoint of 
!    possible separator.
!
!    Input, integer IND, if positive then (XS,YS) has index IND in PVL; else
!    (XS,YS) is on edge joining vertices with indices -IND
!    and SUCC(-IND) in PVL.
!
!    Input, double precision ALPHA, the polar angle of (XS,YS) with respect
!    to (XR,YR).
!
!    Input, double precision THETA, the interior angle at reflex vertex.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer PVL(1:4,1:*), double precision IANG(1:*), the polygon 
!    vertex list, interior angles.
!
!    Output, double precision MINANG, the minimum of the 4 angles in radians.
!
  double precision alpha
  double precision ang
  double precision angle
  double precision beta1
  double precision iang(*)
  integer ind
  integer j
  integer l
  integer, parameter :: loc = 1
  double precision minang
  double precision d_pi
  integer pvl(4,*)
  integer, parameter :: succ = 3
  double precision theta
  double precision vcl(2,*)
  double precision xr
  double precision xs
  double precision yr
  double precision ys
!
  if ( ind > 0 ) then
    j = pvl(succ,ind)
    ang = iang(ind)
  else
    j = pvl(succ,-ind)
    ang = d_pi ( )
  end if

  l = pvl(loc,j)
  beta1 = angle ( xr, yr, xs, ys, vcl(1,l), vcl(2,l) )

  minang = min ( alpha, theta - alpha, ang - beta1, beta1 )

  return
end
subroutine mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )
!
!*******************************************************************************
!
!! MMASEP chooses the best of four separators by the max-min angle criterion.
!
!
!  Purpose: 
!
!    Find best of four possible separators according to
!    max-min angle criterion.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGTOL, the angle tolerance parameter (in 
!    radians) for accepting separator.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the coordinates of 
!    polygon vertices in counter clockwise order where NVRT is number of 
!    vertices; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, integer INDPVL(0:NVRT), indices in PVL of vertices; INDPVL(I) = -K if
!    (XC(I),YC(I)) is extra vertex inserted on edge from K to PVL(SUCC,K).
!
!    Input, double precision IANG(1:*), the interior angle array.
!
!    Input, integer V(1:2), W(1:2), indices in XC, YC in range 0 to NVRT-1; 
!    four possible separators are V(I),W(J), I,J = 1,2.
!
!    Output, integer I1, I2, indices in range 0 to NVRT-1 of best separator
!    according to max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
  double precision alpha
  double precision angle
  double precision angmax
  double precision angmin
  double precision angtol
  double precision beta
  double precision delta
  double precision gamma
  integer i
  integer i1
  integer i2
  double precision iang(*)
  integer indpvl(0:*)
  integer j
  integer k
  integer l
  integer m
  double precision d_pi
  double precision tol
  integer v(2)
  integer w(2)
  double precision xc(0:*)
  double precision yc(0:*)
!
  tol = 100.0D+00 * epsilon ( tol )
  angmax = 0.0D+00

  do i = 1, 2

    l = v(i)
    k = indpvl(l)

    if ( k > 0 ) then
      alpha = iang(k)
    else
      alpha = d_pi ( )
    end if

    do j = 1, 2

      m = w(j)

      if ( l /= m ) then

        k = indpvl(m)

        if ( k > 0 ) then
          beta = iang(k)
        else
          beta = d_pi ( )
        end if

        gamma = angle ( xc(m), yc(m), xc(l), yc(l), xc(l+1), yc(l+1) )

        delta = angle ( xc(l), yc(l), xc(m), yc(m), xc(m+1), yc(m+1) )

        angmin = min ( gamma, alpha-gamma, delta, beta-delta )

        if ( angmin > angmax ) then
          angmax = angmin
          i1 = l
          i2 = m
        end if

      end if

    end do

  end do

  if ( angmax < angtol ) then
    i1 = -1
  end if

  return
end
subroutine mtredg ( utype, i1, i2, i3, ibndry, nt, til, tedg )
!
!*******************************************************************************
!
!! MTREDG sets fields for a triangle as needed by routine TMERGE.
!
!
!  Purpose: 
!
!    Set fields for triangle as needed by routine TMERGE.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical UTYPE, is .TRUE. iff triangle contains two 'U' vertices.
!
!    Input, integer I1, I2, I3, the indices of 3 triangle vertices in VCL; 
!    the first two indices also belong to the next merge edge.
!
!    Input, integer IBNDRY, the index of boundary edge for TEDG.
!
!    Input/output, integer NT, the number of entries in TIL, TEDG so far.
!
!    Input/output, integer TIL(1:NT), the triangle incidence list.
!
!    Input/output, TEDG(1:NT), the triangle edge indices; see routine TMERGE.
!
  integer i1
  integer i2
  integer i3
  integer ibndry
  integer nt
  integer tedg(3,*)
  integer til(3,*)
  logical utype
!
  nt = nt + 1
  til(1,nt) = i1
  til(2,nt) = i2
  til(3,nt) = i3
  tedg(1,nt) = nt

  if ( utype ) then
    tedg(2,nt) = nt - 1
    tedg(3,nt) = ibndry
  else
    tedg(2,nt) = ibndry
    tedg(3,nt) = nt - 1
  end if

  return
end
function prime ( k )
!
!*******************************************************************************
!
!! PRIME returns a prime greater than a given integer K.
!
!
!  Purpose: 
!
!    Return a prime >= K (if possible) from internal array
!    of primes.  More primes can be added if desired.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer K, a positive integer.
!
!    Output, integer PRIME, the smallest prime >= K from internal array 
!    (or largest in array).
!
  integer, parameter :: nprime = 150
!
  integer k
  integer l
  integer m
  integer prime
  integer primes(nprime)
  integer u
!
  save primes

  data primes/17,31,47,61,79,97,113,127,149,163,179,193,211,227,241, &
        257,271,293,307,331,353,379,401,431,457,479,503,541,563,587, &
        613,641,673,701,727,751,773,797,821,853,877,907,929,953,977, &
        1009,1049,1087,1123,1163,1201,1237,1277,1319,1361,1399,1433, &
        1471,1511,1543,1579,1613,1657,1699,1741,1783,1831,1873,1931, &
        1973,2017,2069,2129,2203,2267,2333,2389,2441,2503,2557,2609, &
        2663,2719,2789,2851,2917,2999,3061,3137,3209,3299,3371,3449, &
        3527,3613,3697,3779,3863,3947,4049,4211,4421,4621,4813,5011, &
        5227,5413,5623,5813,6011,6211,6421,6619,6823,7013,7211,7411, &
        7621,7817,8011,8219,8419,8623,8819,9011,9221,9413,9613,9811, &
        10037,10211,10427,10613,10831,11027,11213,11411,11617,11813, &
        12011,12211,12413,12611,12821,13033,13217,13411,13613,13829, &
        14011/

!
  if ( k <= primes(1) ) then
    prime = primes(1)
    return
  else if ( k >= primes(nprime) ) then
    prime = primes(nprime)
    return
  end if
!
!  Use binary search to find prime >= K.
!
  l = 1
  u = nprime

  do

    m = ( l + u ) / 2

    if ( k < primes(m) ) then
      u = m - 1
    else if ( k > primes(m) ) then
      l = m + 1
    else
      prime = primes(m)
      return
    end if

    if ( l > u ) then
      exit
    end if

  end do

  prime = primes(u+1)

  return
end
function d_pi ( )
!
!*******************************************************************************
!
!! D_PI returns the value of pi.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision D_PI, the value of pi.
!
  double precision d_pi
!
  d_pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
subroutine prmdf2 ( ipoly, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, listev )
!
!*******************************************************************************
!
!! PRMDF2 preprocesses a mesh distribution function evaluation.
!
!
!  Purpose: 
!
!    Preprocessing step for evaluating mesh distribution
!    function in polygon IPOLY.  The edges and vertices for
!    which distances must be computed are determined.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer IPOLY, the index of the polygon.
!
!    Input, double precision WSQ, the square of the width of polygon IPOLY.
!
!    Input, integer IVRT(1:*), the indices of polygon vertices in 
!    VCL, ordered by polygon.
!
!    Input, integer XIVRT(1:*), pointers to first vertex of each polygon 
!    in IVRT; vertices of polygon IPOLY are IVRT(I) for I from
!    XIVRT(IPOLY) to XIVRT(IPOLY+1)-1.
!
!    Input, double precision EDGVAL(1:*), a value associated with each edge 
!    of the decomposition.
!
!    Input, double precision VRTVAL(1:*), a value associated with each vertex 
!    of the decomposition.
!
!    Output, integer NEV, the number of edges and vertices for which distances
!    must be evaluated.
!
!    Output, integer IFV, the index of first vertex XIVRT(IPOLY) if LISTEV(NEV)
!    = XIVRT(IPOLY+1) - 1; 0 otherwise.
!
!    Output, integer LISTEV(1:*), an array of length 
!    <= [XIVRT(IPOLY+1)-XIVRT(IPOLY)]*2, containing indices of edges and 
!    vertices mentioned above; indices of vertices are negated.
!
  double precision edgval(*)
  integer i
  integer ifv
  integer im1
  integer ipoly
  integer ivrt(*)
  integer j
  integer l
  integer listev(*)
  integer nev
  double precision vrtval(*)
  double precision wsq
  integer xivrt(*)
!
  ifv = 0
  nev = 0
  im1 = xivrt(ipoly+1) - 1
  l = im1

  do i = xivrt(ipoly), l

    j = ivrt(i)

    if ( vrtval(j) < min ( edgval(i), edgval(im1) ) ) then
      nev = nev + 1
      listev(nev) = -j
    end if

    if ( edgval(i) < wsq ) then
      nev = nev + 1
      listev(nev) = i
    end if

    im1 = i

  end do

  if ( nev > 0 ) then
    if ( listev(nev) == l ) then
      ifv = xivrt(ipoly)
    end if
  end if

  return
end
subroutine ptpolg ( dim, ldv, nv, inc, pgind, vcl, pt, nrml, dtol, inout )
!
!*******************************************************************************
!
!! PTPOLG determines if a point is in, on or outside a polygon.
!
!
!  Purpose: 
!
!    Determine whether a point lies inside, outside, or on
!    boundary of a planar polygon in 2 or 3 dimensional space.
!    It is assumed that point lies in plane of polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer DIM, the dimension of the polygon (2 or 3).
!
!    Input, integer LDV, the leading dimension of VCL array in calling routine.
!
!    Input, integer NV, the number of vertices in polygon.
!
!    Input, integer INC, the increment for PGIND indicating indices of polygon.
!
!    Input, integer PGIND(0:NV*INC), indices in VCL of polygon vertices are in
!    PGIND(0), PGIND(INC), ..., PGIND(NV*INC) with first and
!    last vertices identical.
!
!    Input, double precision VCL(1:DIM,1:*), the vertex coordinate list.
!
!    Input, double precision PT(1:DIM), the point for which in/out test is
!    applied.
!
!    Input, double precision NRML(1:3), the unit normal vector of plane 
!    containing polygon, with vertices oriented counter clockwise with
!    respect to the normal (used iff DIM = 3);
!    The normal is assumed to be (0,0,1) if DIM = 2.
!
!    Input, double precision DTOL, an absolute tolerance to determine 
!    whether a point is on a line or plane.
!
!    Output, integer INOUT, point PT is:
!    +1, inside the polygon, 
!     0, on boundary of polygon, 
!    -1, outside polygon;
!    -2 if error in input parameters.
!
  integer dim
  integer ldv
!
  double precision cp(3)
  double precision de(3)
  double precision dir(3)
  double precision dist
  double precision dotp
  double precision dtol
  integer i
  integer inc
  integer inout
  integer j
  integer k
  integer l
  integer la
  integer lb
  double precision len1
  double precision len2
  integer m
  integer n
  double precision nr(4)
  double precision nrml(3)
  integer nv
  integer pgind(0:*)
  double precision pt(dim)
  double precision rhs(3)
  integer s
  integer sa
  integer sb
  double precision t
  double precision ta
  double precision tol
  double precision vcl(ldv,*)
!
  tol = 100.0D+00 * epsilon ( tol )

  inout = -2

  if ( dim < 2 .or. dim > 3 ) then
    return
  end if
!
!  Direction of ray is from PT through midpoint of first edge
!  such that PT is not collinear with edge. NR is normal of plane
!  containing ray, which is also orthogonal to NRML.
!
  i = 0
  lb = pgind(0)

10 continue

    i = i + 1

    if ( i >= nv ) then
      return
    end if

    la = lb
    lb = pgind(i*inc)

    do j = 1, dim
      de(j) = vcl(j,lb) - vcl(j,la)
      dir(j) = pt(j) - vcl(j,la)
    end do

    if ( dim == 2 ) then
      len1 = de(1)**2 + de(2)**2
      len2 = dir(1)**2 + dir(2)**2
    else
      len1 = de(1)**2 + de(2)**2 + de(3)**2
      len2 = dir(1)**2 + dir(2)**2 + dir(3)**2
    end if

    if ( len1 == 0.0D+00 ) then
      go to 10
    else if ( len2 == 0.0D+00 ) then
      inout = 0
      return
    end if

    if ( dim == 2 ) then
      dotp = abs ( de(1) * dir(1) + de(2) * dir(2)) / sqrt(len1*len2)
    else if ( dim == 3 ) then
      dotp = abs ( de(1) * dir(1) + de(2) * dir(2) + de(3) * dir(3) ) &
        / sqrt(len1*len2)
    end if

    if ( dotp >= 1.0D+00 - tol ) then
      go to 10
    end if

    if ( dim == 2 ) then
      dir(1) = 0.5D+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
      dir(2) = 0.5D+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
      dist = sqrt ( dir(1)**2 + dir(2)**2 )
      dir(1) = dir(1) / dist
      dir(2) = dir(2) / dist
      dir(3) = 0.0D+00
      nr(1) = -dir(2)
      nr(2) = dir(1)
      nr(3) = 0.0D+00
      nr(4) = nr(1) * pt(1) + nr(2) * pt(2)
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
    else if ( dim == 3 ) then
      dir(1) = 0.5D+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
      dir(2) = 0.5D+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
      dir(3) = 0.5D+00 * ( vcl(3,la) + vcl(3,lb) ) - pt(3)
      dist = sqrt ( dir(1)**2 + dir(2)**2 + dir(3)**2 )

      dir(1) = dir(1) / dist
      dir(2) = dir(2) / dist
      dir(3) = dir(3) / dist

      nr(1) = nrml(2)*dir(3) - nrml(3)*dir(2)
      nr(2) = nrml(3)*dir(1) - nrml(1)*dir(3)
      nr(3) = nrml(1)*dir(2) - nrml(2)*dir(1)
      nr(4) = nr(1)*pt(1) + nr(2)*pt(2) + nr(3)*pt(3)
      dist = nr(1)*vcl(1,lb)+nr(2)*vcl(2,lb)+nr(3)*vcl(3,lb) - nr(4)
    end if

    if ( dist > 0.0D+00 ) then
      sb = 1
    else
      sb = -1
    end if

    m = 1
    if ( abs ( dir(2) ) > abs ( dir(1) ) ) then
      m = 2
    end if

    if ( abs ( dir(3) ) > abs ( dir(m) ) ) then
      m = 3
    end if

    k = 1
!
!  For remaining edges of polygon, check whether ray intersects edge.
!  Vertices or edges lying on ray are handled by looking at preceding
!  and succeeding edges not lying on ray.
!
  n = i
  i = i + 1

30 continue

   la = lb
   lb = pgind(i*inc)
   sa = sb

   if ( dim == 2 ) then
     dist = nr(1) * vcl(1,lb) + nr(2)*vcl(2,lb) - nr(4)
   else if ( dim == 3 ) then
     dist = nr(1) * vcl(1,lb) + nr(2)*vcl(2,lb) + nr(3)*vcl(3,lb)- nr(4)
   end if

   if ( abs ( dist) <= dtol ) then
     sb = 0
   else if ( dist > 0.0D+00 ) then
     sb = 1
   else
     sb = -1
   end if

   s = sa * sb

   if ( s < 0 ) then

      if ( dim == 2 ) then
         de(1) = vcl(1,la) - vcl(1,lb)
         de(2) = vcl(2,la) - vcl(2,lb)
         rhs(1) = vcl(1,la) - pt(1)
         rhs(2) = vcl(2,la) - pt(2)
         t = ( rhs(1) * de(2) - rhs(2) * de(1) ) &
           / ( dir(1) * de(2) - dir(2) * de(1) )
      else if ( dim == 3 ) then
         de(1) = vcl(1,la) - vcl(1,lb)
         de(2) = vcl(2,la) - vcl(2,lb)
         de(3) = vcl(3,la) - vcl(3,lb)
         rhs(1) = vcl(1,la) - pt(1)
         rhs(2) = vcl(2,la) - pt(2)
         rhs(3) = vcl(3,la) - pt(3)
         cp(1) = dir(2) * de(3) - dir(3) * de(2)
         cp(2) = dir(3) * de(1) - dir(1) * de(3)
         cp(3) = dir(1) * de(2) - dir(2) * de(1)

         l = 1
         if ( abs ( cp(2) ) > abs ( cp(1) ) ) then
           l = 2
         end if
         if ( abs ( cp(3) ) > abs ( cp(l) ) ) then
           l = 3
         end if

         if ( l == 1 ) then
           t = ( rhs(2) * de(3) - rhs(3) * de(2) ) / cp(1)
         else if ( l == 2 ) then
           t = ( rhs(3) * de(1) - rhs(1) * de(3) ) / cp(2)
         else
           t = ( rhs(1) * de(2) - rhs(2) * de(1) ) / cp(3)
         end if

      end if

      if ( t > dtol ) then
         k = k + 1
      else if ( t >= -dtol ) then
         inout = 0
         return
      end if

   else if ( s == 0 ) then

      l = lb
40    continue
      i = i + 1
      if ( i > nv) i = 1
      if ( i == n) return
      la = lb
      lb = pgind(i*inc)

    if ( dim == 2 ) then
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
    else
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) + nr(3) * vcl(3,lb) - nr(4)
    end if

    if ( abs ( dist) <= dtol ) then
      go to 40
    else if ( dist > 0.0D+00 ) then
      sb = 1
    else
      sb = -1
    end if

    t = ( vcl(m,l) - pt(m) ) / dir(m)

    if ( abs ( t) <= dtol ) then
      inout = 0
      return
    end if

    if ( la /= l ) then
      ta = ( vcl(m,la) - pt(m) ) / dir(m)
      if ( abs ( ta ) <= dtol .or. t * ta < 0.0D+00 ) then
        inout = 0
        return
      end if
    end if

    if ( sa * sb < 0 .and. t > 0.0D+00 ) then
      k = k + 1
    end if

  end if

  i = i + 1

  if ( i > nv ) then
    i = 1
  end if

  if ( i /= n ) then
    go to 30
  end if
!
!  Point lies inside polygon if number of intersections K is odd.
!
  if ( mod ( k, 2 ) == 1 ) then
    inout = 1
  else
    inout = -1
  end if

  return
end
function radians_to_degrees ( angle )
!
!*******************************************************************************
!
!! RADIANS_TO_DEGREES converts an angle from radians to degrees.
!
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision ANGLE, an angle in radians.
!
!    Output, double precision RADIANS_TO_DEGREES, the equivalent angle
!    in degrees.
!
  double precision angle
  double precision d_pi
  double precision radians_to_degrees
!
  radians_to_degrees = ( angle / d_pi ( ) ) * 180.0D+00

  return
end
subroutine randpt ( k, n, seed, axis, nptav, scale, trans, lda, a )
!
!*******************************************************************************
!
!! RANDPT generates N random K-dimensional points from the uniform distribution.
!
!
!  Purpose: 
!
!    Generate N random K-dimensional points from the uniform distribution.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer K, the dimension of points.
!
!    Input, integer N, the number of random points.
!
!    Input/output, integer SEED, a seed for pseudo random number generator.
!
!    Input, integer AXIS, integer NPTAV; if AXIS < 1 or > K, then uniform 
!    random points are generated; if 1 <= AXIS <= K, then an average of NPTAV
!    uniform random points are generated with the same AXIS
!    coordinate on about N/NPTAV random parallel hyperplanes.
!
!    Input, double precision SCALE(K), TRANS(K), the scale and 
!    translation factors for coordinates 1 to K; the I-th coordinate of 
!    random point is R*SCALE(I) + TRANS(I) where 0 < R < 1.
!
!    Input, integer LDA, the leading dimension of array A in calling
!    routine; should be >= K.
!
!    Output, double precision A(LDA,N), an array of N uniform random 
!    K-dimensional points.
!
  integer k
  integer lda
  integer n
!
  double precision a(lda,n)
  integer axis
  integer i
  integer j
  integer m
  integer nptav
  double precision r
  double precision scale(k)
  integer seed
  double precision trans(k)
  real urand
!
  if ( axis < 1 .or. axis > k ) then

    do j = 1, n
      do i = 1, k
        a(i,j) = trans(i) + scale(i) * urand ( seed )
      end do
    end do

  else

    m = int ( urand ( seed ) * 2.0D+00 * nptav + 0.5D+00 )
    r = urand ( seed ) * scale(axis) + trans(axis)

    do j = 1, n

      do i = 1, k
        if ( i == axis ) then
          a(i,j) = r
        else
          a(i,j) = urand ( seed ) * scale(i) + trans(i)
        end if
      end do

      m = m - 1

      if ( m <= 0 ) then
        m = int ( urand ( seed ) * 2.0D+00 * nptav + 0.5D+00 )
        r = urand ( seed ) * scale(axis) + trans(axis)
      end if

    end do

  end if

  return
end
subroutine resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw,  &
  maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )
!
!*******************************************************************************
!
!! RESVRT resolves a reflex vertex of a simply connected polygon.
!
!
!  Purpose: 
!
!    Resolve reflex vertex of simply connected polygon with
!    one or two separators. The reflex vertex must be a 'simple'
!    vertex of the polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer VR, the index in PVL of reflex vertex.
!
!    Input, double precision ANGSPC, the angle spacing parameter used in
!    controlling the vertices to be considered as an endpoint of a separator.
!
!    Input, double precision ANGTOL, the angle tolerance parameter used in 
!    accepting separator(s).
!
!    Input/output, integer NVC, the number of positions used in VCL array.
!
!    Input/output, integer NVERT, the number of positions used in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, integer MAXPV, the maximum size available for PVL array.
!
!    Input, integer MAXIW, the maximum size available for IWK array; should 
!    be about 3 times number of vertices in polygon.
!
!    Input, integer MAXWK, the maximum size available for WK array; should 
!    be about 5 times number of vertices in polygon.
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Output, integer W1, the index in PVL of vertex which is the endpoint 
!    of separator in inner cone or right cone with respect to the reflex vertex.
!
!    Output, integer W2, is 0 if there is only one separator; else index 
!    in PVL of vertex which is endpoint of second separator in left cone.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 5, 6, 7, 206, 207, 208, 209, 210, or 212
!
  integer maxiw
  integer maxpv
  integer maxvc
  integer maxwk
!
  double precision angsep
  double precision angspc
  double precision angtol
  double precision iang(maxpv)
  integer ierror
  integer i
  integer i1
  integer i2
  integer ivis
  integer ivor
  integer ivrt
  integer iwk(maxiw)
  integer l
  integer, parameter :: loc = 1
  integer maxn
  integer nvc
  integer nvert
  integer nvis
  integer nvor
  integer nvrt
  integer nvsvrt
  integer pvl(4,maxpv)
  integer, parameter :: succ = 3
  integer theta
  integer v
  double precision vcl(2,maxvc)
  integer vr
  integer w1
  integer w2
  double precision wk(maxwk)
  integer wkang
  integer xc
  double precision xr
  integer xvor
  integer yc
  double precision yr
  integer yvor
!
!  Determine number of vertices in polygon containing reflex vertex.
!
  ierror = 0
  nvrt = 0
  v = vr

  do

    v = pvl(succ,v)

    if ( v == vr ) then
      exit
    end if

    nvrt = nvrt + 1

  end do

  maxn = nvrt + int ( iang(vr) / angspc )
  l = pvl(loc,vr)
  xr = vcl(1,l)
  yr = vcl(2,l)
!
!  Set up work arrays for routine VISPOL, and determine whether there
!  is enough workspace. XC, YC are d.p. arrays of length NVRT in WK,
!  used for the coordinates of the polygon containing the reflex
!  vertex. MAXN positions are reserved for XC, YC since this is the
!  maximum space required by routine VISVRT. IVIS is an integer array
!  of length MAXN in IWK. IVRT is an integer array of length NVRT in
!  IWK used temporarily for storing indices of vertices in PVL.
!
  if ( maxn + nvrt > maxiw ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESVRT - Fatal error!'
    write ( *, '(a)' ) '  MAXN + NVRT > MAXIW.'
    ierror = 6
    return
  end if

  if ( maxn + maxn > maxwk ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESVRT - Fatal error!'
    write ( *, '(a)' ) '  MAXN + MAXN > MAXWK.'
    ierror = 7
    return
  end if

  ivis = 1
  ivrt = ivis + maxn
  xc = 1
  yc = xc + maxn
  v = pvl(succ,vr)

  do i = 0, nvrt-1
    l = pvl(loc,v)
    wk(xc+i) = vcl(1,l)
    wk(yc+i) = vcl(2,l)
    iwk(ivrt+i) = v
    v = pvl(succ,v)
  end do

  call vispol ( xr, yr, nvrt-1, wk(xc), wk(yc), nvis, iwk(ivis), ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  XC, YC now contain visibility polygon coordinates. Update MAXN
!  and set up d.p. array THETA of length MAXN in WK for routine
!  VISVRT. Elements of IVIS are changed to indices of PVL after call.
!
  maxn = maxn - nvrt + nvis + 1
  theta = yc + maxn

  if ( theta + maxn - 1 > maxwk ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESVRT - Fatal error!'
    write ( *, '(a)' ) '  THETA + MAXN - 1 > MAXWK.'
    ierror = 7
    return
  end if

  call visvrt ( angspc, xr, yr, nvis, wk(xc), wk(yc), iwk(ivis), maxn-1, &
    nvsvrt, wk(theta) )

  wk(theta+nvsvrt) = iang(vr)

  do i = ivis, ivis+nvsvrt
    l = iwk(i)
    if ( l >= 0 ) then
      iwk(i) = iwk(ivrt+l)
    else
      iwk(i) = -iwk(ivrt-l-1)
    end if
  end do
!
!  XC, YC now contain coordinates of visible vertices to be considered
!  as an endpoint of a separator. Set up work arrays for routine
!  VORNBR. Integer array IVOR and d.p. arrays XVOR, YVOR, each of
!  length NVSVRT+1, are added at the end of IWK and WK arrays.
!
  ivor = ivis + nvsvrt + 1
  xvor = theta + nvsvrt + 1
  yvor = xvor + nvsvrt + 1

  if ( ivor + nvsvrt > maxiw ) then
    ierror = 6
    return
  end if

  if ( yvor + nvsvrt > maxwk ) then
    ierror = 7
    return
  end if

  call vornbr ( xr, yr, nvsvrt, wk(xc), wk(yc), nvor, iwk(ivor), wk(xvor), &
    wk(yvor), ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Set up the array WKANG of length NVOR+1 <= NVSVRT+1 in WK for
!  routine FNDSEP. Only Voronoi neighbors are considered as an
!  endpoint of a separator in the first call to FNDSEP. If the
!  minimum angle created at the boundary by the separator(s) is too
!  small, then a second call is made to FNDSEP in which all visible
!  vertices are considered as an endpoint of a separator.
!
  wkang = xvor
  if ( iwk(ivor+nvor) == nvsvrt ) then
    nvor = nvor - 1
  end if

  if ( iwk(ivor) == 0 ) then
    ivor = ivor + 1
    nvor = nvor - 1
  end if

  call fndsep ( angtol+angtol, xr, yr, nvsvrt, wk(xc), wk(yc), iwk(ivis), &
    wk(theta), nvor, iwk(ivor), vcl, pvl, iang, angsep, i1, i2, wk(wkang) )

  if ( angsep < angtol ) then

    ivrt = ivis + nvsvrt + 1

    do i = 1, nvsvrt-1
      iwk(ivrt+i-1) = i
    end do

    call fndsep ( angtol+angtol, xr, yr, nvsvrt, wk(xc), wk(yc), iwk(ivis), &
      wk(theta), nvsvrt-2, iwk(ivrt), vcl, pvl, iang, angsep, i1, i2, wk(wkang) )

  end if
!
!  Insert endpoint(s) of separator(s) in vertex coordinate list and
!  polygon vertex list data structures, if they are not yet there.
!
  if ( i2 == -1 ) then
    w2 = 0
  else if ( iwk(ivis+i2) < 0 ) then
    call insvr2 ( wk(xc+i2), wk(yc+i2), -iwk(ivis+i2), nvc, nvert, maxvc, &
      maxpv, vcl, pvl, iang, w2, ierror )
    if ( ierror /= 0 ) then
      return
    end if
  else
    w2 = iwk(ivis+i2)
  end if

  if ( iwk(ivis+i1) < 0 ) then
    call insvr2 ( wk(xc+i1), wk(yc+i1), -iwk(ivis+i1), nvc, nvert, maxvc, &
      maxpv, vcl, pvl, iang, w1, ierror )
    if ( ierror /= 0 ) then
      return
    end if
  else
    w1 = iwk(ivis+i1)
  end if

  return
end
subroutine rotiar ( n, arr, shift )
!
!*******************************************************************************
!
!! ROTIAR rotates the elements of an integer array.
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer N, the number of elements of array.
!
!    Input/output, integer ARR(0:N-1), the array to be shifted.
!
!    Input, integer SHIFT, the amount of (left) shift or rotation; 
!    ARR(SHIFT) on input becomes ARR(0) on output.
!
  integer n
!
  integer a
  integer arr(0:n-1)
  integer b
  integer i
  integer j
  integer k
  integer l
  integer m
  integer r
  integer sh
  integer shift
  integer t
!
  sh = mod ( shift, n )

  if ( sh < 0 ) then
    sh = sh + n
  end if

  if ( sh == 0 ) then
    return
  end if

  a = n
  b = sh

  do
    r = mod ( a, b )
    a = b
    b = r
    if ( r <= 0 ) then
      exit
    end if
  end do

  m = n / a - 1

  do i = 0, a-1

    t = arr(i)
    k = i

    do j = 1, m
      l = k + sh
      if ( l >= n ) then
        l = l - n
      end if
      arr(k) = arr(l)
      k = l
    end do

    arr(k) = t

  end do

  return
end
subroutine rotipg ( xeye, yeye, nvrt, xc, yc, ierror )
!
!*******************************************************************************
!
!! ROTIPG rotates the vertex indices of a simple polygon.
!
!
!  Purpose: 
!
!    Rotate the indices of the vertices of a simple polygon
!    and possibly insert one vertex so that (XC(0),YC(0)) is the
!    point on the horizontal line through (XEYE,YEYE) and on the
!    boundary of the polygon which is closest to and to the right
!    of (XEYE,YEYE). (XEYE,YEYE) is an eyepoint in the interior or
!    blocked exterior of the polygon. In the former (latter) case,
!    the vertices must be in counter clockwise (CW) order.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XEYE, YEYE, the coordinates of eyepoint.
!
!    Input/output, integer NVRT, the number of vertices on boundary of simple 
!    polygon.  On output, NVRT is increased by 1 if the closest vertex
!    is a new vertex.
!
!    Input/output, double precision XC(0:NVRT), YC(0:NVRT), the vertices 
!    of polygon in counter clockwise (or clockwise) order if eyepoint is 
!    interior (or blocked exterior);
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT))
!    On output, the polygon vertices in same orientation
!    as input but with indices rotated and possibly
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)) has been added.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 205.
!
  integer nvrt
!
  integer a
  integer b
  double precision dy
  integer i
  integer ierror
  integer irgt
  integer j
  integer k
  integer l
  integer m
  integer n
  integer r
  double precision tol
  double precision xc(0:nvrt+1)
  double precision xeye
  double precision xint
  double precision xrgt
  double precision xt
  double precision yc(0:nvrt+1)
  double precision yeye
  double precision yeyemt
  double precision yeyept
  double precision yt
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  dy = 0.0D+00
  do i = 0, nvrt-1
    dy = max ( dy, abs ( yc(i+1) - yc(i) ) )
  end do

  yeyemt = yeye - tol * dy
  yeyept = yeye + tol * dy
  n = nvrt + 1
  irgt = n
  xrgt = 0.0D+00
!
!  Determine closest point on boundary which is to the right of
!  (XEYE,YEYE) and on the horizontal line through (XEYE,YEYE).
!  The closest point must be on an edge which intersects the
!  horizontal line and has (XEYE,YEYE) to its left.
!
  do i = 0, nvrt-1

    if ( yc(i) > yeyept .or. yc(i+1) < yeyemt ) then
      cycle
    end if

    if ( yc(i) < yeyemt .and. yc(i+1) > yeyept ) then
      xint = ( yeye - yc(i) ) * ( xc(i+1) - xc(i) ) / ( yc(i+1) - yc(i) ) + xc(i)
      if ( xint > xeye ) then
         if ( xint < xrgt .or. irgt == n ) then
            irgt = -(i + 1)
            xrgt = xint
         end if
      end if
    else if ( yc(i) >= yeyemt .and. yc(i+1) > yeyept ) then
      if ( xc(i) > xeye ) then
         if ( xc(i) < xrgt .or. irgt == n ) then
            irgt = i
            xrgt = xc(i)
         end if
      end if
    else if ( yc(i) < yeyemt .and. yc(i+1) <= yeyept ) then
      if ( xc(i+1) > xeye ) then
         if ( xc(i+1) < xrgt .or. irgt == n ) then
            irgt = i + 1
            xrgt = xc(i+1)
         end if
      end if
    end if

  end do

  if ( irgt == n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTIPG - Fatal error!'
    write ( *, '(a)' ) '  IRGT == N.'
    ierror = 205
    return
  end if

  if ( irgt == 0 .or. irgt == nvrt ) then
    return
  end if

  if ( irgt < 0 ) then
    irgt = -irgt
    do i = nvrt, irgt, -1
      xc(i+1) = xc(i)
      yc(i+1) = yc(i)
    end do
    xc(irgt) = xrgt
    yc(irgt) = yeye
    nvrt = nvrt + 1
  end if
!
!  Rotate the indices of the vertices so that (XC(IRGT),YC(IRGT))
!  becomes (XC(0),YC(0)). Compute A = GCD(NVRT,IRGT).
!
  a = nvrt
  b = irgt

  do

    r = mod ( a, b )
    a = b
    b = r

    if ( r <= 0 ) then
      exit
    end if

  end do

  m = nvrt / a - 1

  do i = 0, a-1

    xt = xc(i)
    yt = yc(i)
    k = i

    do j = 1, m
      l = k + irgt
      if ( l >= nvrt ) then
        l = l - nvrt
      end if
      xc(k) = xc(l)
      yc(k) = yc(l)
      k = l
    end do

    xc(k) = xt
    yc(k) = yt

  end do

  xc(nvrt) = xc(0)
  yc(nvrt) = yc(0)

  return
end
subroutine rotpg ( nvrt, xc, yc, i1, i2, ibot, costh, sinth )
!
!*******************************************************************************
!
!! ROTPG rotates a convex polygon.
!
!
!  Purpose: 
!
!    Rotate convex polygon so that a line segment joining two
!    of its vertices is parallel to y-axis.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices on the boundary of
!    the convex polygon.
!
!    Input/output, double precision XC(0:NVRT), YC(0:NVRT).  The vertex 
!    coordinates in counter clockwise order;
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!    On output, the rotated vertex coordinates; indices are
!    also rotated so that (XC(0),YC(0)) = (XC(NVRT),YC(NVRT))
!    is top vertex and (XC(IBOT),YC(IBOT)) is bottom vertex.
!
!    Input, integer I1, I2, the index of vertices of line segment; I1, I2 > 0.
!
!    Output, integer IBOT, the index of bottom vertex.
!
!    Output, double precision COSTH, SINTH, the values COS(THETA) and 
!    SIN(THETA) where THETA in [-PI,PI] is the rotation angle.
!
  integer nvrt
!
  integer a
  integer b
  double precision costh
  integer i
  integer i1
  integer i2
  integer ibot
  integer itop
  integer j
  integer k
  integer l
  integer m
  double precision d_pi
  integer r
  double precision sinth
  double precision theta
  double precision tol
  double precision x0
  double precision xc(0:nvrt)
  double precision y0
  double precision yc(0:nvrt)
!
  tol = 100.0D+00 * epsilon ( tol )

  itop = i1
  ibot = i2

  if ( yc(i1) == yc(i2) ) then
    if ( xc(i1) < xc(i2) ) then
      theta = -d_pi ( ) / 2.0D+00
    else
      theta = d_pi ( ) / 2.0D+00
    end if
  else
    if ( yc(i1) < yc(i2) ) then
      itop = i2
      ibot = i1
    end if
    theta = d_pi ( ) / 2.0D+00 - atan2 ( yc(itop) - yc(ibot), xc(itop) - xc(ibot) )
  end if

  costh = cos(theta)
  sinth = sin(theta)

  do i = 1, nvrt
    x0 = xc(i)
    xc(i) = costh * x0 - sinth * yc(i)
    yc(i) = sinth * x0 + costh * yc(i)
  end do
!
!  Rotate indices.
!
  if ( itop /= nvrt ) then

    a = nvrt
    b = itop

    do

      r = mod ( a, b )
      a = b
      b = r

      if ( r <= 0 ) then
        exit
      end if

    end do

    m = nvrt / a - 1

    do i = 1, a

      x0 = xc(i)
      y0 = yc(i)
      k = i

      do j = 1, m
        l = k + itop
        if ( l > nvrt ) then
          l = l - nvrt
        end if
        xc(k) = xc(l)
        yc(k) = yc(l)
        k = l
      end do

      xc(k) = x0
      yc(k) = y0

    end do

    ibot = ibot - itop
    if ( ibot < 0 ) then
      ibot = ibot + nvrt
    end if

  end if

  xc(0) = xc(nvrt)
  yc(0) = yc(nvrt)

  return
end
subroutine sepmdf ( angtol, nvrt, xc, yc, arpoly, mean, mdftr, indpvl, &
  iang, i1, i2 )
!
!*******************************************************************************
!
!! SEPMDF splits a polygon according to the mesh distribution function.
!
!
!  Purpose: 
!
!    Determine separator to split convex polygon into two
!    parts based on mesh distribution function.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGTOL, the angle tolerance parameter (in radians).
!
!    Input, integer NVRT, the number of vertices in polygon.
!
!    Input, double precision XC(0:NVRT),YC(0:NVRT), the coordinates of polygon
!    vertices in counter clockwise order, translated so that centroid is at 
!    origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, double precision ARPOLY, the area of polygon.
!
!    Input, double precision MEAN, the mean mdf value in polygon.
!
!    Input, double precision MDFTR(0:NVRT-1), the mean mdf value in each 
!    triangle of polygon; triangles are determined by polygon vertices 
!    and centroid.
!
!    Input, integer INDPVL(0:NVRT), the indices in PVL of vertices; 
!    INDPVL(I) = -K if (XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, double precision IANG(1:*), the interior angle array.
!
!    Output, integer I1, I2, indices in range 0 to NVRT-1 of best separator
!    according to MDF and max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
  integer nvrt
!
  double precision angle
  double precision angtol
  double precision areatr
  double precision arpoly
  integer hi
  integer i
  integer i1
  integer i2
  double precision iang(*)
  integer indpvl(0:nvrt)
  integer l
  integer m
  double precision mdftr(0:nvrt-1)
  double precision mean
  double precision d_pi
  double precision sum2
  double precision tol
  integer v(2)
  integer w(2)
  double precision xc(0:nvrt)
  double precision yc(0:nvrt)
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine triangle with highest mean mesh density; then determine
!  triangles adjacent to this triangle with mesh density >= MEAN
!  such that the area of these triangles is <= ARPOLY/2.
!  Note that twice the triangle area is computed.
!
  hi = 0
  do i = 1, nvrt-1
    if ( mdftr(i) > mdftr(hi) ) then
      hi = i
    end if
  end do

  sum2 = xc(hi) * yc(hi+1) - xc(hi+1) * yc(hi)

  l = hi - 1
  if ( l < 0 ) then
    l = nvrt - 1
  end if

  m = hi + 1
  if ( m >= nvrt ) then
    m = 0
  end if

  do

    if ( mdftr(l) >= mdftr(m) ) then
      i = l
    else
      i = m
    end if

    if ( mdftr(i) < mean ) then
      exit
    end if

    areatr = xc(i) * yc(i+1) - xc(i+1) * yc(i)
    sum2 = sum2 + areatr

    if ( sum2 > arpoly ) then
      exit
    end if

    if ( i == l ) then
      l = l - 1
      if ( l < 0 ) then
        l = nvrt - 1
      end if
    else
      m = m + 1
      if ( m >= nvrt ) then
        m = 0
      end if
    end if

  end do

  l = l + 1

  if ( l >= nvrt ) then
    l = 0
  end if
!
!  Interchange role of L and M depending on angle determined by
!  (XC(M),YC(M)), (0,0), and (XC(L),YC(L)).
!  Possible separators are L,M; L,M+1; L+1,M; L+1,M+1.
!
  if ( angle ( xc(m), yc(m), 0.0D+00, 0.0D+00, xc(l), yc(l) ) > d_pi ( ) ) then
    i = l
    l = m
    m = i
  end if

  v(1) = l
  v(2) = l - 1
  if ( v(2) < 0 ) then 
    v(2) = nvrt - 1
  end if

  w(1) = m
  w(2) = m + 1
  if ( w(2) >= nvrt ) then
    w(2) = 0
  end if

  call mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

  return
end
subroutine sepshp ( angtol, nvrt, xc, yc, indpvl, iang, i1, i2, wk, ierror )
!
!*******************************************************************************
!
!! SEPSHP splits a convex polygon according to shape.
!
!
!  Purpose: 
!
!    Determine separator to split convex polygon into two
!    parts based on shape (diameter) of polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGTOL, the angle tolerance parameter (in radians).
!
!    Input, integer NVRT, the number of vertices in polygon.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the coordinates of 
!    polygon vertices in counter clockwise order, translated so that 
!    centroid is at origin;
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, integer INDPVL(0:NVRT), the indices in PVL of vertices;
!    INDPVL(I) = -K if (XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, double precision IANG(1:*), the interior angle array.
!
!    Output, integer I1, I2, the indices in range 0 to NVRT-1 of best separator
!    according to shape and max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
!    Workspace, double precision WK(1:2*NVRT).
!
!    Output, integer IERROR, error flag.  For abnormal return,
!    IERROR is set to 200.
!
  integer nvrt
!
  double precision angtol
  double precision dist
  double precision dx
  double precision dy
  integer i
  integer i1
  integer i2
  double precision iang(*)
  integer ierror
  integer indpvl(0:nvrt)
  integer k
  integer n
  double precision d_pi
  double precision pimtol
  double precision tol
  integer v(2)
  integer w(2)
  double precision wk(2*nvrt)
  double precision xa
  double precision xc(0:nvrt)
  double precision ya
  double precision yc(0:nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine diameter of polygon. Possible separators endpoints (two
!  on each side of polygon) are nearest to perpendicular bisector of
!  diameter. (XA,YA) and (XA+DX,YA+DY) are on bisector. Distance to
!  bisector is proportional to two times triangle area.
!
  pimtol = d_pi ( ) - tol
  n = 0
  do i = 0, nvrt-1
    k = indpvl(i)
    if ( k > 0 ) then
      if ( iang(k) < pimtol ) then
         n = n + 1
         wk(n) = xc(i)
         wk(n+nvrt) = yc(i)
      end if
    end if
  end do

  call diam2 ( n, wk, wk(nvrt+1), i1, i2, dist, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  if ( i1 > i2 ) then
    i = i1
    i1 = i2
    i2 = i
  end if

  dx = wk(i2+nvrt) - wk(i1+nvrt)
  dy = wk(i1) - wk(i2)
  xa = 0.5D+00 * ( wk(i1) + wk(i2) - dx )
  ya = 0.5D+00 * ( wk(i1+nvrt) + wk(i2+nvrt) - dy )

  i = i1 - 1

20 continue

  if ( xc(i) == wk(i1) .and. yc(i) == wk(i1+nvrt) ) then
    i1 = i
  else
    i = i + 1
    go to 20
  end if

  i = max ( i2-1, i1+1 )

30 continue

  if ( xc(i) == wk(i2) .and. yc(i)  ==  wk(i2+nvrt) ) then
    i2 = i
  else
    i = i + 1
    go to 30
  end if

  i = i1 + 1

40 continue

  dist = dx * ( yc(i) - ya ) - dy * ( xc(i) - xa )

  if ( dist >= 0.0D+00 ) then
    v(1) = i - 1
    v(2) = i
  else
    i = i + 1
    go to 40
  end if

  i = i2 + 1

50 continue

  if ( i >= nvrt) then
    i = 0
  end if

  dist = dx * ( yc(i) - ya ) - dy * ( xc(i) - xa )

  if ( dist <= 0.0D+00 ) then
    w(1) = i - 1
    w(2) = i
    if ( i <= 0 ) then
      w(1) = nvrt - 1
    end if
  else
    i = i + 1
    go to 50
  end if

  call mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

  return
end
subroutine sfdwmf ( l, r, psi, indp, loch )
!
!*******************************************************************************
!
!! SFDWMF sifts PSI(INDP(L)) down a heap.
!
!
!  Purpose: 
!
!    Sift PSI(INDP(L)) down heap which has maximum PSI value
!    at root of heap and is maintained by pointers in INDP.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer L, the element of heap to be sifted down.
!
!    Input, integer R, the upper bound of heap.
!
!    Input, double precision PSI(1:*), the key values for heap.
!
!    Input/output, integer INDP(1:R), the indices of PSI which are 
!    maintained in heap.
!
!    Input/output, integer LOCH(1:*), the location of indices in heap 
!    (inverse of INDP).
!
  integer r
!
  integer i
  integer indp(r)
  integer j
  integer k
  integer l
  integer loch(*)
  double precision psi(*)
  double precision t
!
  i = l
  j = 2 * i
  k = indp(i)
  t = psi(k)

  do while ( j <= r )

    if ( j < r ) then
      if ( psi(indp(j)) < psi(indp(j+1)) ) then
        j = j + 1
      end if
    end if

    if ( t >= psi(indp(j)) ) then
      exit
    end if

    indp(i) = indp(j)
    loch(indp(i)) = i
    i = j
    j = 2 * i

  end do

  indp(i) = k
  loch(k) = i

  return
end
subroutine sfupmf ( r, psi, indp, loch )
!
!*******************************************************************************
!
!! SFUPMF sifts PSI(INDP(R)) up a heap.
!
!
!  Purpose: 
!
!    Sift PSI(INDP(R)) up heap which has maximum PSI value
!    at root of heap and is maintained by pointers in INDP.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer R, the element of heap to be sifted up.
!
!    Input, double precision PSI(1:*), the key values for heap.
!
!    Input/output, integer INDP(1:R), the indices of PSI which are 
!    maintained in heap.
!
!    Input/output, integer LOCH(1:*), the location of indices in heap 
!    (inverse of INDP).
!
  integer r
!
  integer i
  integer indp(r)
  integer j
  integer k
  integer loch(*)
  double precision psi(*)
  double precision t
!
  i = r
  j = int ( i / 2 )
  k = indp(i)
  t = psi(k)

  do

    if ( i <= 1 ) then
      exit
    end if

    if ( t <= psi(indp(j)) ) then
      exit
    end if

    indp(i) = indp(j)
    loch(indp(i)) = i
    i = j
    j = int ( i / 2 )

  end do

  indp(i) = k
  loch(k) = i

  return
end
subroutine shrnk2 ( nvrt, xc, yc, sdist, nshr, xs, ys, iedge, ierror )
!
!*******************************************************************************
!
!! SHRNK2 shrinks a convex polygon.
!
!
!  Purpose: 
!
!    Shrink a convex polygon, with vertices given in counter clockwise
!    order and with all interior angles < PI, by distance SDIST(I)
!    for Ith edge, I = 0,...,NVRT-1.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices on the boundary of convex 
!    polygon.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    in counter clockwise order;
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, double precision SDIST(0:NVRT-1), the nonnegative shrink 
!    distances for edges.
!
!    Output, integer NSHR, the number of vertices on boundary of shrunken 
!    polygon; 0 if shrunken polygon is empty else 3 <= NSHR <= NVRT.
!
!    Output, double precision XS(0:NSHR), YS(0:NSHR), the coordinates of
!    shrunken polygon in counter clockwise order if NSHR > 0; 
!    (XS(0),YS(0)) = (XS(NSHR),YS(NSHR)).
!
!    Output, integer IEDGE(0:NVRT), the indices of edges of shrunken polygon in
!    range from 0 to NVRT-1.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 202.
!
  integer nvrt
!
  double precision alpha
  logical first
  integer i
  integer iedge(0:nvrt)
  integer ierror
  integer j
  integer k
  integer lr
  integer lrline
  integer nshr
  logical parall
  double precision d_pi
  double precision pi2
  double precision sdist(0:nvrt-1)
  double precision theta
  double precision tol
  double precision xc(0:nvrt)
  double precision xs(0:nvrt)
  double precision yc(0:nvrt)
  double precision ys(0:nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  pi2 = 2.0D+00 * d_pi ( )
  alpha = atan2 ( yc(1)-yc(0), xc(1)-xc(0) )

  call xline ( xc(0), yc(0), xc(1), yc(1), xc(1), yc(1), xc(2), yc(2), &
    sdist(0), sdist(1), xs(1), ys(1), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  iedge(0) = 0
  iedge(1) = 1
  i = 2
  j = 0
  nshr = 1
  first = .true.
!
!  First while loop processes edges subtending angle <= PI
!  with respect to first edge.
!
10 continue

  theta = atan2 ( yc(i+1)-yc(i), xc(i+1)-xc(i) ) - alpha
  if ( theta < 0.0D+00 ) then
    theta = theta + pi2
  end if

  if ( theta > d_pi ( ) + tol ) then
    go to 40
  end if

20 continue

  lr = lrline ( xs(nshr), ys(nshr), xc(i), yc(i), xc(i+1), yc(i+1), &
    sdist(i) )

  if ( lr < 0 ) then
    go to 30
  end if

  nshr = nshr - 1
  if ( nshr >= 1 ) then
    go to 20
  end if

30 continue

  if ( nshr < 1 .and. abs ( theta - d_pi ( ) ) <= tol ) then
    go to 90
  end if

  k = iedge(nshr)
  nshr = nshr + 1

  call xline ( xc(k), yc(k), xc(k+1), yc(k+1), xc(i), yc(i), xc(i+1), &
    yc(i+1), sdist(k), sdist(i), xs(nshr), ys(nshr), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  iedge(nshr) = i
  i = i + 1
  go to 10
!
!  Second while loop processes remaining edges.
!
40 continue

  if ( first ) then
    first = .false.
    go to 50
  end if

  lr = lrline ( xs(j), ys(j), xc(i), yc(i), xc(i+1), yc(i+1), sdist(i) )

  if ( lr <= 0 ) then
    go to 70
  end if

50 continue

  if ( nshr <= j ) then
    go to 90
  end if

  lr = lrline ( xs(nshr), ys(nshr), xc(i), yc(i), xc(i+1), yc(i+1), sdist(i) )

  if ( lr >= 0 ) then
    nshr = nshr - 1
    go to 50
  end if

  k = iedge(nshr)
  nshr = nshr + 1

  call xline ( xc(k), yc(k), xc(k+1), yc(k+1), xc(i), yc(i), xc(i+1), &
    yc(i+1), sdist(k), sdist(i), xs(nshr), ys(nshr), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  iedge(nshr) = i

60 continue

  lr = lrline ( xs(j+1), ys(j+1), xc(i), yc(i), xc(i+1), yc(i+1), sdist(i) )

  if ( lr >= 0 ) then
    j = j + 1
    go to 60
  end if

  k = iedge(j)
  call xline ( xc(k), yc(k), xc(k+1), yc(k+1), xc(i), yc(i), xc(i+1), &
    yc(i+1), sdist(k), sdist(i), xs(j), ys(j), parall )

  if ( parall ) then
    ierror = 202
    go to 90
  end if

  xs(nshr+1) = xs(j)
  ys(nshr+1) = ys(j)
  iedge(nshr+1) = iedge(j)

70 continue

  i = i + 1

  if ( i < nvrt ) then
    go to 40
  end if

  if ( j > 0 ) then
    do i = 0, nshr+1-j
      xs(i) = xs(i+j)
      ys(i) = ys(i+j)
      iedge(i) = iedge(i+j)
    end do
  end if

  nshr = nshr + 1 - j
  return

90 continue

  nshr = 0

  return
end
subroutine spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc,  &
  maxhv, maxpv, maxiw, maxwk, holv, vcl, regnum, hvl, pvl, iang, iwk, &
  wk, ierror )
!
!*******************************************************************************
!
!! SPDEC2 decomposes a polygonal region with holes into simple polygons.
!
!
!  Discussion: 
!
!    Decompose general polygonal region with interfaces and
!    holes into simple polygons using vertex coordinate list,
!    head vertex list, and polygon vertex list data structures.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGSPC, the angle spacing parameter in radians 
!    used in controlling vertices to be considered as an endpoint of a separator.
!
!    Input, double precision ANGTOL, the angle tolerance parameter in radians 
!    used in accepting separator(s).
!
!    Input/output, integer NVC, the number of vertex coordinates or positions 
!    used in VCL array.
!
!    Input/output, integer NPOLG, the number of polygonal subregions or 
!    positions used in HVL array.
!
!    Input/output, integer NVERT, the number of polygon vertices or positions 
!    used in PVL array.
!
!    Input, integer NHOLE, the number of holes and hole interfaces.
!
!    Input, integer NHOLA, the number of 'attached' holes; these holes are 
!    attached to the outer boundary of a subregion through vertices
!    or cut interfaces and have their edges in consecutive order on the boundary.
!
!    Input, integer MAXVC, the maximum size available for VCL array, should 
!    be >= number of vertex coordinates required for decomposition.
!
!    Input, integer MAXHV, the maximum size available for HVL, REGNUM arrays, 
!    should be >= number of polygons required for decomposition.
!
!    Input, integer MAXPV, the maximum size available for PVL, IANG arrays; 
!    should be >= number of polygon vertices required for decomposition.
!
!    Input, integer MAXIW, the maximum size available for IWK array; should be 
!    about 3 times maximum number of vertices in any polygon.
!
!    Input, integer MAXWK, the maximum size available for WK array; should be 
!    about 5 times maximum number of vertices in any polygon.
!
!    Input, integer HOLV(1:NHOLE*2+NHOLA), the indices in PVL of bottom or top 
!    vertex of holes; first (next) NHOLE entries are for top (bottom)
!    vertices of holes and hole interfaces, with top (bottom)
!    vertices sorted in decreasing (increasing) lexicographic
!    (y,x) order of coord; last NHOLA entries are for attached
!    holes; if bottom vertex of attached hole is a simple
!    vertex of boundary curve containing the hole then entry
!    contains index of bottom vertex otherwise entry contains
!    index of top vertex (which is simple).
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT), 
!    the polygon vertex list and interior angles; see routine DSPGDC for more 
!    details.  Note: The data structures should be as output from routines
!    DSMCPR or DSPGDC.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 4, 5, 6, 7, 206 to 210, 212, 218, or 219.
!
  integer maxhv
  integer maxiw
  integer maxpv
  integer maxvc
  integer maxwk
!
  double precision angspc
  double precision angtol
  logical ci
  logical cj
  integer, parameter :: edgv = 4
  integer holv(*)
  integer hvl(maxhv)
  integer i
  double precision iang(maxpv)
  integer ierror
  integer iwk(maxiw)
  integer j
  integer nhola
  integer nhole
  integer npolg
  integer nvc
  integer nvert
  integer p
  double precision d_pi
  double precision piptol
  integer, parameter :: polg = 2
  integer pvl(4,maxpv)
  integer regnum(maxhv)
  integer, parameter :: succ = 3
  double precision tol
  double precision vcl(2,maxvc)
  integer vr
  integer w1
  integer w2
  double precision wk(maxwk)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  For each simple hole, find cut edge from top vertex of hole to
!  a point on the outer boundary above top vertex, and update
!  VCL, HVL, PVL, IANG.
!
  piptol = d_pi ( ) + tol

  do i = 1, nhole

    call jnhole ( holv(i), angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
      maxwk, vcl, hvl, pvl, iang, iwk, wk, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPDEC2 - Fatal error!'
      write ( *, '(a)' ) '  JNHOLE returned an error condition.'
      return
    end if

  end do
!
!  Resolve remaining vertices in HOLV array if they are reflex
!  vertices. These vertices may no longer be reflex if they are the
!  endpoint of a cut edge from the top vertex of another hole or
!  of a previous separator.
!
  do i = nhole+1, nhole+nhole+nhola

    vr = holv(i)

    if ( iang(vr) > piptol ) then

      call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
        maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call insed2 ( vr, w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
        pvl, iang, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( w2 > 0 ) then
        call insed2 ( vr, w2, npolg, nvert, maxhv, maxpv, &
          vcl, regnum, hvl, pvl, iang, ierror )
        if ( ierror /= 0 ) then
          return
        end if
      end if

    end if

  end do

  if ( nhola == 0 ) then
    return
  end if
!
!  Check that polygons are simple. If polygon is simply-connected and
!  not simple then find a simple reflex vertex in polygon to resolve.
!
  p = 1

30 continue

  if ( p > npolg ) then
    return
  end if

  i = hvl(p)

  do

    if ( pvl(polg,pvl(edgv,i)) == p ) then
      go to 50
    end if

    i = pvl(succ,i)

    if ( i == hvl(p) ) then
      exit
    end if

  end do

  p = p + 1
  go to 30

50 continue

  ci = .true.

  do

    j = pvl(succ,i)
    cj = ( pvl(polg,pvl(edgv,j)) == p )

    if ( .not. ci .and. .not. cj .and. iang(j) > piptol ) then
      exit
    end if

    i = j
    ci = cj

  end do

  vr = j
  call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
    maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call insed2 ( vr, w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
    pvl, iang, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  if ( w2 > 0 ) then

    call insed2 ( vr, w2, npolg, nvert, maxhv, maxpv, &
      vcl, regnum, hvl, pvl, iang, ierror )

    if ( ierror /= 0 ) then
      return
    end if

  end if

  go to 30

end
subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierror )
!
!*******************************************************************************
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2-D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to triangulation.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer I, the index in VCL of the new vertex.
!
!    Input/output, integer TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input, integer MAXST, the maximum size available for the STACK array.
!
!    Input/output, integer BTRI, BEDG; on input, if positive, are the
!    triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, double precision VCL(2,*), the coordinates of the vertices.
!
!    Input/output, integer TIL(3,*), the triangle incidence list.  May be updated
!    on output because of swaps.
!
!    Input/output, integer TNBR(3,*), the triangle neighbor list; negative 
!    values are used for links of the counter-clockwise linked list of boundary 
!    edges; May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(1:MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERROR is set to 8 for abnormal return.
!
  integer maxst
!
  integer a
  integer b
  integer bedg
  integer btri
  integer c
  integer diaedg
  integer e
  integer ee
  integer em1
  integer ep1
  integer f
  integer fm1
  integer fp1
  integer i
  integer ierror
  integer l
  integer r
  integer s
  integer stack(maxst)
  integer swap
  integer t
  integer til(3,*)
  integer tnbr(3,*)
  integer top
  integer tt
  integer u
  double precision vcl(2,*)
  double precision x
  double precision y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  ierror = 0
  x = vcl(1,i)
  y = vcl(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( til(1,t) == i ) then
      e = 2
      b = til(3,t)
    else if ( til(2,t) == i ) then
      e = 3
      b = til(1,t)
    else
      e = 1
      b = til(2,t)
    end if

    a = til(e,t)
    u = tnbr(e,t)

    if ( tnbr(1,u) == t ) then
      f = 1
      c = til(3,u)
    else if ( tnbr(2,u) == t ) then
      f = 2
      c = til(1,u)
    else
      f = 3
      c = til(2,u)
    end if

    swap = diaedg ( x, y, vcl(1,a), vcl(2,a), vcl(1,c), vcl(2,c), &
      vcl(1,b), vcl(2,b) )

    if ( swap == 1 ) then

      em1 = i_wrap ( e - 1, 1, 3 )
      ep1 = i_wrap ( e + 1, 1, 3 )
      fm1 = i_wrap ( f - 1, 1, 3 )
      fp1 = i_wrap ( f + 1, 1, 3 )
      
      til(ep1,t) = c
      til(fp1,u) = i
      r = tnbr(ep1,t)
      s = tnbr(fp1,u)
      tnbr(ep1,t) = u
      tnbr(fp1,u) = t
      tnbr(e,t) = s
      tnbr(f,u) = r

      if ( tnbr(fm1,u) > 0 ) then
        top = top + 1
        stack(top) = u
      end if

      if ( s > 0 ) then

        if ( tnbr(1,s) == u ) then
          tnbr(1,s) = t
        else if ( tnbr(2,s) == u ) then
          tnbr(2,s) = t
        else
          tnbr(3,s) = t
        end if

        top = top + 1

        if ( top > maxst ) then
          ierror = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( tnbr(ee,tt) > 0 )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == a ) then
            ee = 3
          else if ( til(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

      if ( r > 0 ) then

        if ( tnbr(1,r) == t ) then
          tnbr(1,r) = u
        else if ( tnbr(2,r) == t ) then
          tnbr(2,r) = u
        else
          tnbr(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( tnbr(ee,tt) > 0 )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == b ) then
            ee = 3
          else if ( til(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

    end if

  end do

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tmerge ( inter, nbl, ncr, chbl, chcr, ldv, vcl, til, tedg, &
  ierror )
!
!*******************************************************************************
!
!! TMERGE forms triangles near the boundary by merging vertex chains.
!
!
!  Purpose: 
!
!    Form triangles in strip near boundary of polygon or
!    inside polygon by merging two chains of vertices.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical INTER, is .TRUE. iff at least one interior mesh vertex.
!
!    Input, integer NBL, the number of vertices on boundary cycle if INTER,
!    otherwise on left boundary chain.
!
!    Input, integer NCR, the number of vertices on closed walk if INTER,
!    otherwise on right boundary chain.
!
!    Input, integer CHBL(0:NBL), the indices in VCL of vertices on boundary cycle
!    or left boundary chain; if INTER, CHBL(NBL) = CHBL(0).
!
!    Input, integer CHCR(0:NCR), the indices in VCL of vertices on closed walk
!    or right boundary chain; if INTER, CHCR(NCR) = CHCR(0),
!    otherwise CHCR(0) is not referenced.
!
!    Input, integer LDV, the leading dimension of VCL in calling routine.
!
!    Input, double precision VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, integer TIL(1:3,1:NT), the triangle incidence list, where NT =
!    NBL + NCR - K where K = 0 if INTER, else K = 2.
!
!    Output, integer TEDG(1:3,1:NT), the TEDG(J,I) refers to edge with vertices
!    TIL(J:J+1,I) and contains index of merge edge or NBL+NCR+1 for edge of 
!    chains.  Note: It is assumed there is enough space in 2 arrays.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 230.
!
  integer ldv
  integer nbl
  integer ncr
!
  integer chbl(0:nbl)
  integer chcr(0:ncr)
  integer diaedg
  integer i
  integer ibndry
  integer ierror
  integer in
  logical inter
  integer j
  integer lri
  integer lrip1
  integer lrline
  integer nl
  integer nr
  integer nt
  integer tedg(3,nbl+ncr)
  integer til(3,nbl+ncr)
  double precision vcl(ldv,*)
  double precision xi
  double precision xip1
  double precision xj
  double precision xjp1
  double precision yi
  double precision yip1
  double precision yj
  double precision yjp1
!
  ierror = 0
  ibndry = nbl + ncr + 1
  nt = 0

  if ( inter ) then

    nl = nbl
    nr = ncr
    i = 0
    j = 0

  else

    call mtredg ( .true., chbl(1), chcr(1), chbl(0), ibndry, nt, til, tedg )

    tedg(2,1) = ibndry
    if ( nbl + ncr <= 3 ) then
      return
    end if

    nl = nbl - 1
    nr = ncr - 1
    i = 1
    j = 1
    lri = 1
    lrip1 = 1

  end if
!
!  Main while loop for determining next triangle and edge.
!
10 continue

  if ( i >= nl .or. j >= nr ) then
    go to 20
  end if

   xi = vcl(1,chbl(i))
   yi = vcl(2,chbl(i))
   xip1 = vcl(1,chbl(i+1))
   yip1 = vcl(2,chbl(i+1))
   xj = vcl(1,chcr(j))
   yj = vcl(2,chcr(j))
   xjp1 = vcl(1,chcr(j+1))
   yjp1 = vcl(2,chcr(j+1))
   in = diaedg ( xjp1, yjp1, xj, yj, xi, yi, xip1, yip1 )

   if ( inter ) then
     lri = lrline ( xi, yi, xj, yj, xjp1, yjp1, 0.0D+00 )
     lrip1 = lrline ( xip1, yip1, xj, yj, xjp1, yjp1, 0.0D+00 )
   end if

   if ( in <= 0 .or. lri <= 0 .and. lrip1 <= 0 ) then

     call mtredg ( .true., chbl(i+1), chcr(j), chbl(i), ibndry, nt, til, tedg )
     i = i + 1

   else

     call mtredg ( .false., chbl(i), chcr(j+1), chcr(j), ibndry, nt, til, tedg )
     j = j + 1

   end if

  go to 10
!
!  Add remaining triangles at end of strip or bottom of polygon.
!
20 continue

  if ( i < nl ) then

    if ( .not. inter .and. j == nr ) then
      nl = nl + 1
    end if

    do

      call mtredg ( .true., chbl(i+1), chcr(j), chbl(i), ibndry, nt, til, tedg )
      i = i + 1

      if ( i >= nl ) then
        exit
      end if

    end do

  else
!
!  J < NR .OR. I = NL = J = NR = 1
!
    if ( .not. inter .and. i == nl ) then
      nr = nr + 1
    end if

40  continue

    call mtredg ( .false., chbl(i), chcr(j+1), chcr(j), ibndry, nt, til, tedg )

    if ( inter ) then

      lri = lrline ( vcl(1,chbl(i)), vcl(2,chbl(i)), &
        vcl(1,chcr(j+1)), vcl(2,chcr(j+1)), vcl(1,chcr(j)), &
        vcl(2,chcr(j)), 0.0D+00 )

      if ( lri >= 0 ) then
        ierror = 230
        return
      end if

    end if

    j = j + 1

    if ( j < nr ) then
      go to 40
    end if

  end if

  if ( inter ) then
    if ( tedg(2,1) == 0 ) then
      tedg(2,1) = nbl + ncr
    else
      tedg(3,1) = nbl + ncr
    end if
  end if

  return
end
subroutine trinbr ( nvc, ntri, til, tnbr, htsiz, maxedg, ht, edge, ierror )
!
!*******************************************************************************
!
!! TRINBR determines the neighboring triangles of every triangle.
!
!
!  Purpose: 
!
!    Determine the neighboring triangle, if any, along each edge
!    of every triangle of triangulation of polygonal region.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVC, the number of vertices in triangulation.
!
!    Input, integer NTRI, the number of triangles in triangulation.
!
!    Input, integer TIL(1:3,1:NTRI), the triangle incidence list; TIL(1:3,I) 
!    contains indices in VCL of 3 vertices of Ith triangle in counter 
!    clockwise order.
!
!    Input, integer HTSIZ, the size of hash table HT; should be a prime number
!    which is about NB where NB is number of boundary edges.
!
!    Input, integer MAXEDG, the maximum size available for EDGE array; should
!    be about 2*NB.
!
!    Output, integer TNBR(1:3,1:NTRI), the triangle neighbor list; positive 
!    elements are indices of TIL; zero elements indicate boundary edges.
!
!    Workspace, integer HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG), the hash table 
!    and edge records used to determine matching occurrences of triangle edges
!    by calling routine EDGHT.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 1.
!
  integer htsiz
  integer maxedg
  integer ntri
!
  integer e
  integer edge(4,maxedg)
  integer hdfree
  integer ht(0:htsiz-1)
  integer i
  integer i3
  integer ierror
  integer j
  integer, save, dimension ( 3 ) :: jp1 = (/ 2, 3, 1 /)
  integer last
  integer nvc
  integer t
  integer til(3,ntri)
  integer tnbr(3,ntri)
  integer w
!
  ierror = 0
  hdfree = 0
  last = 0
  ht(0:htsiz-1) = 0

  do i = 1, ntri

    i3 = i * 3

    do j = 1, 3

      call edght ( til(j,i), til(jp1(j),i), i3+j-1, nvc, htsiz, maxedg, &
        hdfree, last, ht, edge, w, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( w > 0 ) then
        t = w / 3
        e = mod ( w, 3 ) + 1
        tnbr(e,t) = i
        tnbr(j,i) = t
      end if

    end do

  end do

  do while ( hdfree /= 0 )
    edge(1,hdfree) = 0
    hdfree = edge(4,hdfree)
  end do

  do i = 1, last

    if ( edge(1,i) /= 0 ) then
      t = edge(3,i)/3
      e = mod ( edge(3,i), 3 ) + 1
      tnbr(e,t) = 0
    end if

  end do

  return
end
subroutine tripr2 ( nvc, npolg, nvert, maxvc, maxti, maxiw, maxwk, h, vcl,  &
  hvl, pvl, iang, ntri, til, vstart, vnum, tstart, iwk, wk, ierror )
!
!*******************************************************************************
!
!! TRIPR2 generates triangles inside each convex polygon of a decomposition.
!
!
!  Purpose: 
!
!    Generate mesh vertices and triangles inside each convex
!    polygon of decomposition according to mesh spacings in H array
!    to get a triangulation of a polygonal region.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input/output, integer NVC, the number of vertex coordinates or positions
!    used in VCL array.
!
!    Input, integer NPOLG, the number of polygonal subregions or positions 
!    used in HVL array.
!
!    Input, integer NVERT, the number of polygon vertices or positions used 
!    in PVL array.
!
!    Input, integer MAXVC, the maximum size available for VCL array, should 
!    be >= number of mesh vertices in triangulation of region.
!
!    Input, integer MAXTI, the maximum size available for TIL array, should 
!    be >= number of triangles in triangulation of region.
!
!    Input, integer MAXIW, the maximum size available for IWK array, should 
!    be >= 5*(NBC+NCW)+2 where NBC is maximum number of mesh edges
!    on boundary of a polygon, NCW is maximum number of edges
!    on boundary of interior triangulation.
!
!    Input, integer MAXWK, the maximum size available for WK array, should 
!    be >= 5*NVRT+4 where NVRT is max no. of vertices in a polygon.
!
!    Input, double precision H(1:NPOLG), the mesh spacings for the polygons 
!    of the decomposition.
!
!    Input/output, double precision VCL(2,MAXVC), the vertex coordinates.
!
!    Input, integer HVL(1:NPOLG), the head vertex list.
!
!    Input, integer PVL(1:4,1:NVERT), double precision IANG(1:NVERT), the 
!    polygon vertex list and interior angles; see routine DSPGDC for 
!    more details.
!
!    Output, integer NTRI, the number of triangles in triangulation of region.
!
!    Output, integer TIL(1:3,1:NTRI), the triangle incidence list; TIL(1:3,I) 
!    contains indices in VCL of 3 vertices of Ith triangle in counter
!    clockwise order.
!
!    Output, integer VSTART(1:NVERT), the start location in VCL for mesh 
!    vertices on each edge in PVL if there are any, else 0.
!
!    Output, integer VNUM(1:NVERT), the number of mesh vertices on interior
!    of each edge in PVL; entry is negated if mesh vertices are
!    listed in backward order in VCL.
!
!    Output, integer TSTART(1:NPOLG), the start location in TIL of triangles in
!    each polygon; TIL(1:3,I) for I=TSTRT(K),...,TSTRT(K+1)-1
!    are the triangles in the K-th polygon.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 6, 7, 9, 10, 200, 202, 230, or 231
!
  integer maxiw
  integer maxti
  integer maxvc
  integer maxwk
  integer npolg
  integer nvert
!
  integer bndcyc
  double precision h(npolg)
  integer hvl(npolg)
  integer i
  double precision iang(nvert)
  integer ierror
  integer iwk(maxiw)
  integer j
  integer k
  integer, parameter :: loc = 1
  integer nbc
  integer ntri
  integer nvc
  integer nvrt
  double precision d_pi
  double precision pimtol
  integer pvl(4,nvert)
  integer, parameter :: succ = 3
  integer til(3,maxti)
  double precision tol
  integer tstart(npolg)
  double precision vcl(2,maxvc)
  integer vnum(nvert)
  integer vstart(nvert)
  double precision wk(maxwk)
  integer xc
  integer yc
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  ntri = 0
  pimtol = d_pi ( ) - tol


  call bedgmv ( nvc, npolg, nvert, maxvc, h, vcl, hvl, pvl, vstart, vnum, &
    ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPR2 - Fatal error!'
    write ( *, '(a)' ) '  Error return from BEDGMV.'
    return
  end if

  do k = 1, npolg

    nvrt = 0
    nbc = 0
    i = hvl(k)

    do

      if ( iang(i) < pimtol ) then
        nvrt = nvrt + 1
      endif

      nbc = nbc + 1 + abs ( vnum(i))
      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    if ( nbc + 1 > maxiw ) then
      ierror = 6
      return
    end if

    if ( 2*nvrt + 2 > maxwk ) then
      ierror = 7
      return
    end if

    xc = 1
    yc = xc + nvrt + 1
    bndcyc = 1

    do

      j = pvl(loc,i)

      if ( iang(i) < pimtol ) then
        wk(xc) = vcl(1,j)
        wk(yc) = vcl(2,j)
        xc = xc + 1
        yc = yc + 1
      end if

      iwk(bndcyc) = j
      bndcyc = bndcyc + 1

      if ( vnum(i) >= 0 ) then
        do j = vstart(i), vstart(i)+vnum(i)-1
          iwk(bndcyc) = j
          bndcyc = bndcyc + 1
        end do
      else
        do j = vstart(i)-vnum(i)-1, vstart(i),-1
          iwk(bndcyc) = j
          bndcyc = bndcyc + 1
        end do
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    wk(xc) = wk(1)
    wk(yc) = wk(nvrt+2)
    iwk(bndcyc) = iwk(1)
    xc = 1
    yc = xc + nvrt + 1
    bndcyc = 1
    tstart(k) = ntri + 1
!
!  Generate a Delaunay triangulation inside a convex polygon.
!
    call trpolg ( nvrt, wk(xc), wk(yc), h(k), nbc, iwk(bndcyc), 2, nvc, ntri, &
      maxvc, maxti, maxiw-nbc-1, maxwk-2*nvrt-2, vcl, til, iwk(nbc+2), &
      wk(2*nvrt+3), ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIPR2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from TRPOLG.'
      return
    end if

  end do

  return
end
subroutine trisiz ( ntrid, npolg, hvl, pvl, area, psi, h, indp, loch )
!
!*******************************************************************************
!
!! TRISIZ smooths the mean mesh distribution function.
!
!
!  Purpose: 
!
!    Smooth PSI (mean mesh distribution function) values using
!    heap so that they differ by a factor of at most 4 in adjacent
!    polygons and then compute triangle sizes for each polygon.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NTRID, the desired number of triangles in mesh.
!
!    Input, integer NPOLG, the number of polygons or positions used in HVL array.
!
!    Input, integer HVL(1:NPOLG), the head vertex list.
!
!    Input, integer PVL(1:4,1:*), the polygon vertex list.
!
!    Input, double precision AREA(1:NPOLG), the area of convex polygons 
!    in decomposition.
!
!    Input/output, double precision PSI(1:NPOLG), the mean mdf values in 
!    the convex polygons.
!
!    Output, double precision H(1:NPOLG), the triangle size for convex polygons.
!
!    Workspace, integer INDP(1:NPOLG), the indices of polygon or PSI which 
!    are maintained in heap according to PSI values.
!
!    Workspace, integer LOCH(1:NPOLG), the location of polygon indices in heap.
!
  integer npolg
!
  double precision area(npolg)
  integer, parameter :: edgv = 4
  double precision factor
  double precision h(npolg)
  integer hvl(npolg)
  integer i
  integer indp(npolg)
  integer j
  integer k
  integer l
  integer loch(npolg)
  integer ntrid
  integer, parameter :: polg = 2
  double precision psi(npolg)
  integer pvl(4,*)
  integer r
  integer, parameter :: succ = 3
!
  factor = 0.25D+00

  call ivec_identity ( npolg, indp )
  call ivec_identity ( npolg, loch )

  k = int ( npolg / 2 )

  do l = k, 1, -1
    call sfdwmf ( l, npolg, psi, indp, loch )
  end do

  do r = npolg, 2, -1

    j = indp(1)
    indp(1) = indp(r)
    loch(indp(1)) = 1
    call sfdwmf ( 1, r-1, psi, indp, loch )
    i = hvl(j)

    do

      k = pvl(edgv,i)

      if ( k > 0 ) then

        k = pvl(polg,k)

        if ( psi(k) < psi(j) * factor ) then
          psi(k) = psi(j) * factor
          call sfupmf ( loch(k), psi, indp, loch )
        end if

      end if

      i = pvl(succ,i)

      if ( i == hvl(j) ) then
        exit
      end if

    end do

  end do

  psi(1:npolg) = psi(1:npolg) / dot_product ( psi(1:npolg), area(1:npolg) )

  h(1:npolg) = sqrt ( 2.0D+00  / ( dble ( ntrid ) * psi(1:npolg) ) )

  return
end
subroutine trpolg ( nvrt, xc, yc, h, nbc, bndcyc, ldv, nvc, ntri, maxvc,  &
  maxti, maxiw, maxwk, vcl, til, iwk, wk, ierror )
!
!*******************************************************************************
!
!! TRPOLG generates a Delaunay triangular mesh inside a convex polygon.
!
!
!  Discussion:
!
!    A quasi-uniform grid of spacing H is used.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the vertex coordinates 
!    in counter clockwise order; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)); it is 
!    assumed that all interior angles are < PI.
!
!    Input, double precision H, the spacing of mesh vertices in polygon.
!
!    Input, integer NBC, the size of BNDCYC.
!
!    Input/output, integer BNDCYC(0:NBC), the indices in VCL of mesh 
!    vertices of boundary cycle; BNDCYC(0) = BNDCYC(NBC); contains (XC(I),YC(I)).
!
!    Input, integer LDV, the leading dimension of VCL in calling routine.
!
!    Input/output, integer NVC, the number of coordinates or positions used 
!    in VCL array.
!
!    Input/output, integer NTRI, the number of triangles or positions used 
!    in TIL.
!
!    Input, integer MAXVC, the maximum size available for VCL array.
!
!    Input, integer MAXTI, the maximum size available for TIL array.
!
!    Input, integer MAXIW, the maximum size available for IWK array, should 
!    be >= 6*(1 + INT(DIAM/H)) + 4*(NBC + NCW) where DIAM is
!    diameter of polygon, NCW is number of edges on boundary
!    of interior triangulation.
!
!    Input, integer MAXWK, the maximum size available for WK array, should 
!    be >= 3*NVRT+2.
!
!    Input/output, double precision VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer TIL(1:3,1:NTRI), the triangle incidence list.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, double precision WK(1:MAXWK).
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 3, 6, 7, 9, 10, 200, 202, 230, or 231.
!
  integer ldv
  integer maxiw
  integer maxti
  integer maxvc
  integer maxwk
  integer nbc
  integer nvrt
!
  integer bndcyc(0:nbc)
  double precision costh
  integer cwalk
  double precision dist
  double precision h
  double precision hs
  integer i
  integer i1
  integer i2
  integer ibot
  integer iedge
  integer ierror
  integer ind
  logical inter
  integer iwk(maxiw)
  integer maxcw
  integer mbc
  integer ncw
  integer nshr
  integer nt
  integer ntri
  integer nvc
  integer sdist
  double precision sinth
  double precision smdist
  integer sptr
  integer tedg
  integer til(3,maxti)
  double precision vcl(ldv,maxvc)
  double precision wk(maxwk)
  double precision x0
  double precision xc(0:nvrt)
  double precision xi
  integer xs
  double precision y0
  double precision yc(0:nvrt)
  double precision yi
  double precision yr
  integer ys
!
  ierror = 0

  if ( nvrt + 1 > maxiw ) then
    ierror = 6
    return
  end if

  if ( 3*nvrt + 2 > maxwk ) then
    ierror = 7
    return
  end if

  xs = 1
  ys = xs + nvrt + 1
  sdist = ys + nvrt + 1
  iedge = 1
  hs = h / sqrt ( 2.0D+00 )
  wk(sdist:sdist+nvrt-1) = hs

  call shrnk2 ( nvrt, xc, yc, wk(sdist), nshr, wk(xs), wk(ys), iwk(iedge), &
    ierror )

  if ( ierror /= 0 ) then
    return
  end if

  inter = ( nshr > 0 )

  if ( inter ) then

    call diam2 ( nshr, wk(xs+1), wk(ys+1), i1, i2, dist, ierror )

    if ( ierror /= 0 ) then
      return
    end if

    call rotpg ( nshr, wk(xs), wk(ys), i1, i2, ibot, costh, sinth )

    maxcw = 6 * ( 1 + int ( ( wk(ys) - wk(ys+ibot) ) / h ) )

    if ( maxcw + 1 > maxiw ) then
      ierror = 6
      return
    end if

    cwalk = 1

    call inttri ( nshr, wk(xs), wk(ys), h, ibot, costh, sinth, ldv, nvc, ntri, &
      maxvc, maxti, maxcw, vcl, til, ncw, iwk(cwalk), ierror )

    if ( ierror /= 0 ) then
      return
    end if
!
!  Determine the mesh vertex which should be moved to front of
!  BNDCYC - closest to CWALK(0) and also with y-coordinate >
!  that of CWALK(0) when rotated if NCW > 0.
!
    x0 = vcl(1,iwk(cwalk))
    y0 = vcl(2,iwk(cwalk))

    if ( ncw > 0 ) then
      yr = sinth * x0 + costh * y0
    end if

    smdist = 100000.0D+00 * h**2

    do i = 0, nbc-1

      xi = vcl(1,bndcyc(i))
      yi = vcl(2,bndcyc(i))

      if ( ncw > 0 ) then
        if ( sinth * xi + costh * yi <= yr ) then
          cycle
        end if
      end if

      dist = ( xi - x0 )**2 + ( yi - y0 )**2

      if ( dist < smdist ) then
        smdist = dist
        ind = i
      end if

    end do

    call rotiar ( nbc, bndcyc, ind )
    bndcyc(nbc) = bndcyc(0)
    nt = nbc + ncw
    tedg = cwalk + ncw + 1

  else

    call diam2 ( nvrt, xc(1), yc(1), i1, i2, dist, ierror )

    if ( ierror /= 0 ) then
      return
    end if

    ind = 0

    do

      if ( ind >= nbc ) then
        exit
      end if

      if ( xc(i1) == vcl(1,bndcyc(ind)) .and. &
           yc(i1) == vcl(2,bndcyc(ind)) ) then
        exit
      end if

      ind = ind + 1

    end do

    call rotiar ( nbc, bndcyc, ind )
    bndcyc(nbc) = bndcyc(0)
    mbc = 1

    do

      if ( mbc >= nbc ) then
        exit
      end if

      if ( xc(i2) == vcl(1,bndcyc(mbc)) .and. &
           yc(i2) == vcl(2,bndcyc(mbc)) ) then
        exit
      end if

      mbc = mbc + 1

    end do

    ind = nbc

    do i = mbc+1, mbc+(nbc-mbc-1)/2
      ind = ind - 1
      i1 = bndcyc(i)
      bndcyc(i) = bndcyc(ind)
      bndcyc(ind) = i1
    end do

    bndcyc(nbc) = bndcyc(mbc)
    nt = nbc - 2
    tedg = 1
!
!  Left boundary chain contains mesh vertices BNDCYC(0:MBC)
!  and right chain contains BNDCYC(0,MBC+1:NBC); MBC < NBC.
!
  end if

  if ( ntri + nt > maxti ) then
    ierror = 9
    return
  else if ( tedg + 4*nt - 1 > maxiw ) then
    ierror = 6
    return
  end if

  if ( inter ) then
    call tmerge ( inter, nbc, ncw, bndcyc, iwk(cwalk), ldv, vcl, &
      til(1,ntri+1), iwk(tedg), ierror )
  else
    call tmerge ( inter, mbc, nbc-mbc, bndcyc, bndcyc(mbc), ldv, vcl, &
      til(1,ntri+1), iwk(tedg), ierror )
  end if

  if ( ierror /= 0 ) then
    return
  end if

  sptr = tedg + 3 * nt

  call cvdtri ( inter, ldv, nt, vcl, til(1,ntri+1), iwk(tedg), iwk(sptr), &
    ierror )

  ntri = ntri + nt

  return
end
function umdf2 ( x, y )
!
!*******************************************************************************
!
!! UMDF2 is a dummy mesh distribution function.
!
!
!  Purpose: 
!
!    Dummy user-supplied mesh distribution function which
!    is provided if heuristic mesh distribution function is used.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a point.
!
!    Output, double precision UMDF2, the mesh distribution function value 
!    at (X,Y)
!
  double precision umdf2
  double precision x
  double precision y
!
  umdf2 = 1.0D+00

  return
end
function urand ( iy )
!
!*******************************************************************************
!
!! URAND is a uniform random number generator.
!
!
!  Discussion:
!
!    URAND is a uniform random number generator based on theory and
!    suggestions given in D. E. Knuth (1969), Vol. 2. The integer IY
!    should be initialized to an arbitrary integer prior to the first
!    call to URAND. The calling program should not alter the value of
!    IY between subsequent calls to URAND. Values of URAND will be
!    returned in the interval (0,1).
!
!  Reference:
!
!    Forsythe, Malcolm, Moler, 
!    page 246.
!
!  Modified:
!
!    12 July 1999
!
!  Parameters:
!
!    Input, integer IY, the seed value.
!
!    Output, double precision URAND, the random value.
!
  double precision halfm
  integer, save :: ia = 0
  integer, save :: ic = 0
  integer, parameter :: itwo = 2
  integer iy
  integer m
  integer, save :: m2 = 0
  integer, save :: mic = 0
  real, save :: s = 0.0E+00
  real urand
!
!  If first entry, compute machine integer word length.
!
  if ( m2 == 0 ) then

    m = 1

    do

      m2 = m
      m = itwo * m2

      if ( m <= m2 ) then
        exit
      end if

    end do

    halfm = m2
!
!  Compute multiplier and increment for linear congruential method.
!
    ia = 8 * int ( halfm * atan ( 1.0D+00 ) / 8.0D+00 ) + 5
    ic = 2 * int ( halfm * ( 0.5D+00 - sqrt ( 3.0D+00 ) / 6.0D+00 ) ) + 1
    mic = ( m2 - ic ) + m2
!
!  S is the scale factor for converting to floating point.
!
    s = 0.5D+00 / halfm

  end if
!
!  Compute next random number.
!
  iy = iy * ia
!
!  The following statement is for computers which do not allow
!  integer overflow on addition.
!
  if ( iy > mic ) then
    iy = (iy - m2) - m2
  end if

  iy = iy + ic
!
!  The following statement is for computers where the word
!  length for addition is greater than for multiplication.
!
  if ( iy / 2 > m2 ) then
    iy = (iy - m2) - m2
  end if
!
!  The following statement is for computers where integer
!  overflow affects the sign bit.
!
  if ( iy < 0 ) then
    iy = ( iy + m2 ) + m2
  end if

  urand = real ( iy ) * s

  return
end
subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )
!
!*******************************************************************************
!
!! VBEDG determines visible boundary edges of a 2D triangulation.
!
!
!  Purpose: 
!
!    Determine boundary edges of 2-D triangulation which are
!    visible from point (X,Y) outside convex hull.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a 2-D point outside
!    the convex hull.
!
!    Input, double precision VCL(1:2,1:*), the coordinates of 2-D vertices.
!
!    Input, integer TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer TNBR(1:3,1:*), the triangle neighbor list; negative 
!    values are used for links of counter clockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer LTRI, LEDG.  On input, if LTRI /= 0 then they 
!    are assumed to be as defined below and are not changed, else they are 
!    updated.  On output, LTRI is the index of the boundary triangle to the
!    left of leftmost boundary triangle visible from (X,Y), and LEDG is the
!    boundary edge of triangle LTRI to left of leftmost
!    boundary edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer RTRI, on input, the index of boundary triangle 
!    to begin search at.  On output, the index of rightmost boundary triangle 
!    visible from (X,Y).
!
!    Input/output, integer REDG.  On input, the edge of triangle RTRI that 
!    is visible from (X,Y).  On output, REDG has been updated so that this
!    is still true. 1 <= REDG <= 3.
!
  integer a
  integer b
  integer e
  integer i_wrap
  integer l
  logical ldone
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer redg
  integer rtri
  integer t
  integer til(3,*)
  integer tnbr(3,*)
  double precision vcl(2,*)
  double precision x
  double precision y
!
!  Find rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

10 continue

  l = -tnbr(redg,rtri)
  t = l / 3
  e = mod ( l, 3 ) + 1
  a = til(e,t)

  if ( e <= 2 ) then
    b = til(e+1,t)
  else
    b = til(1,t)
  end if

  lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), 0.0D+00 )

  if ( lr > 0 ) then
    rtri = t
    redg = e
    go to 10
  end if

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = til(e,t)
    e = i_wrap ( e-1, 1, 3 )

    do while ( tnbr(e,t) > 0 )

      t = tnbr(e,t)

      if ( til(1,t) == b ) then
        e = 3
      else if ( til(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = til(e,t)
    lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
subroutine vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis, ierror )
!
!*******************************************************************************
!
!! VISPOL computes the visibility polygon.
!
!
!  Purpose: 
!
!    Compute the visibility polygon VP from an eyepoint in
!    the interior or blocked exterior of a simple polygon P or
!    on the boundary of a simply connected polygonal region P.
!    In the latter case, the interior angles at all vertices must
!    be strictly between 0 and 2*PI.
!
!  Discussion:
!
!    On input, XC and YC contain vertex coordinates of P. During
!    the algorithm, part of XC, YC is used as a stack, which, on
!    output, contains the vertex coordinates of VP. The stack
!    vertices overwrite the input vertices as the input vertices
!    are scanned. Elements of IVIS are set when vertices are added
!    to the stack; these values may have +NV or -NV added to them
!    to indicate that stack point has same angle as previous one.
!
!  Reference:
!
!    Barry Joe and R. B. Simpson, 
!    BIT 27 (1987), pages 458-473.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XEYE, YEYE, the coordinates of eyepoint; must 
!    be a simple vertex if it lies on the boundary (i.e. occurs only once).
!
!    Input, integer NVRT, the upper subscript of XC, YC (approximate
!    number of vertices).
!
!    Input/output, double precision XC(0:NVRT), YC(0:NVRT).  On input, if 
!    eyepoint is interior or blocked exterior then arrays contain coordinates 
!    in counter clockwise or clockwise order, respectively, with 
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)); (XC(0),YC(0)) is a vertex visible from
!    (XEYE,YEYE), e.g. as computed by routine ROTIPG.  If eyepoint is a vertex
!    of P then arrays contain coordinates in counter clockwise order; 
!    (XC(0),YC(0)) is successor vertex of (XEYE,YEYE); (XC(NVRT),YC(NVRT)) is
!    predecessor vertex of (XEYE,YEYE).
!    On output, XC and YC contain the vertices of VP in counter clockwise order;
!    if eyepoint is interior or blocked exterior then
!    (XC(0),YC(0)) = (XC(NVIS),YC(NVIS)), else (XC(0),YC(0))
!    and (XC(NVIS),YC(NVIS)) are the successor and
!    predecessor vertices of (XEYE,YEYE) in VP.
!
!    Output, integer NVIS, the upper subscript of XC, YC on output (approximate
!    number of vertices of VP); NVIS <= NVRT.
!
!    Output, integer IVIS(0:NVIS), contains information about the vertices 
!    of VP with respect to the vertices of P; IVIS(I) = K if (XC(I),YC(I))
!    is the vertex of index K in the input polygon; IVIS(I) = -K if 
!    (XC(I),YC(I)) is on the interior of the edge joining vertices of index 
!    K-1 and K in input polygon
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 206, 207, 208, 209, or 210
!
  integer nvrt
!
  logical beye
  integer cur
  integer i
  integer ierror
  integer ivis(0:nvrt)
  integer lr
  integer lrline
  integer nv
  integer nvis
  integer oper
  integer top
  double precision xc(0:nvrt)
  double precision xe
  double precision xeye
  double precision xw
  double precision yc(0:nvrt)
  double precision ye
  double precision yeye
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!     Variables in common block GVPVAR:
!        NV - NVRT
!        OPER - operation code 1 to 7 for LEFT, RIGHT, SCANA, SCANB,
!              SCANC, SCAND, FINISH
!        CUR - index of current vertex of P in XC, YC arrays
!        TOP - index of top vertex of stack in XC, YC arrays
!              (TOP <= CUR is always satisfied)
!        XE,YE - XEYE,YEYE
!        XW,YW - coordinates of point on last or second-last edge
!              processed (needed for routines VPSCNB, VPSCNC, VPSCND)
!        BEYE - .TRUE. iff eyepoint is on boundary
!
  ierror = 0
  beye = xc(0) /= xc(nvrt) .or. yc(0) /= yc(nvrt)
  nv = nvrt
  xe = xeye
  ye = yeye
  ivis(0) = 0
  cur = 1

  if ( beye ) then

    do

      lr = lrline ( xc(nv-1), yc(nv-1), xe, ye, xc(nv), yc(nv), 0.0D+00 )

      if ( lr /= 0 ) then
        exit
      end if
      nv = nv - 1

    end do

  end if

  do

    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(0), yc(0), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    cur = cur + 1

  end do

  if ( lr == -1 ) then
    oper = 1
    if ( cur == 1 ) then
      top = 1
      ivis(1) = cur
    else if ( beye ) then
      top = 1
      xc(0) = xc(cur-1)
      yc(0) = yc(cur-1)
      ivis(0) = cur - 1
      xc(1) = xc(cur)
      yc(1) = yc(cur)
      ivis(1) = cur
    else
      top = 2
      xc(1) = xc(cur-1)
      yc(1) = yc(cur-1)
      ivis(1) = cur - 1 + nv
      xc(2) = xc(cur)
      yc(2) = yc(cur)
      ivis(2) = cur
    end if
  else
    oper = 3
    top = 0
    if ( beye .and. cur > 1 ) then
      xc(0) = xc(cur-1)
      yc(0) = yc(cur-1)
      ivis(0) = cur - 1
    end if
  end if
!
!  Angular displacement of stack points are in nondecreasing order,
!  with at most two consecutive points having the same displacement.
!
30 continue

   if ( oper == 1 ) then
     call vpleft ( xc, yc, ivis )
   else if ( oper == 2 ) then
     call vprght ( xc, yc, ivis, ierror )
   else if ( oper == 3 ) then
     call vpscna ( xc, yc, ivis, ierror )
   else if ( oper == 4 ) then
     call vpscnb ( xc, yc, ivis, ierror )
   else if ( oper == 5 ) then
     call vpscnc ( xc, yc, ivis, ierror )
   else
     call vpscnd ( xc, yc, ivis, ierror )
   end if

   if ( ierror /= 0 ) then
     nvis = top
     return
   end if

  if ( oper <= 6 ) then
    go to 30
  end if
!
!  Add or subtract NV from those IVIS values which are used to
!  indicate that stack point has same angle as previous one.
!
  do i = 1, top

    if ( ivis(i) > nv ) then
      ivis(i) = ivis(i) - nv
    else if ( ivis(i) < -nv ) then
      ivis(i) = ivis(i) + nv
    end if

  end do

  nvis = top

  return
end
subroutine visvrt ( angspc, xeye, yeye, nvis, xc, yc, ivis, maxn, nvsvrt, theta )
!
!*******************************************************************************
!
!! VISVRT determines a list of visible vertices.
!
!
!  Purpose: 
!
!    Determine a list of visible vertices, ordered by
!    increasing "polar angle", on the boundary of the visibilty
!    polygon from boundary eyepoint (XEYE,YEYE).  This list
!    includes the vertices of visibility polygon such that a
!    line segment from (XEYE,YEYE) to vertex lies in interior
!    of polygon, as well as extra points on edges which subtend
!    an angle >= 2*ANGSPC at (XEYE,YEYE).  These extra points are
!    at an equal angular spacing of >= ANGSPC and < 2*ANGSPC. The
!    successor and predecessor of eyepoint are included in list.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision ANGSPC, the angle spacing parameter in radians 
!    which controls how many extra points become visible vertices.
!
!    Input, double precision XEYE, YEYE, the coordinates of boundary eyepoint.
!
!    Input, integer NVIS, (number of vertices of visibility polygon) - 2.
!
!    Input/output, double precision XC(0:NVIS), YC(0:NVIS), on input, the 
!    coordinates of the vertices of visibility polygon in counter clockwise 
!    order; (XC(0),YC(0)) and (XC(NVIS),YC(NVIS)) are the successor and 
!    predecessor vertices of eyepoint in visibility polygon; at most 2
!    consecutive vertices have same polar angle with respect to eyepoint.
!    On output, coordinates of visible vertices which overwrite the 
!    input coordinates.
!
!    Input/output, IVIS(0:NVIS), on input, contains information about the 
!    vertices of XC, YC arrays with respect to the original polygon from
!    which visibility polygon is computed; if IVIS(I) >= 0
!    then (XC(I),YC(I)) has index I in original polygon;
!    if IVIS(I) < 0 then (XC(I),YC(I)) is on the edge
!    ending at vertex of index -IVIS(I) in original polygon;
!    indexing starts at 0 from successor of eyepoint.
!    On output, coordinates of visible vertices
!    which overwrite the input coordinates.
!
!    Input, integer MAXN, the upper bound on NVSVRT; should be at least
!    NVIS + INT(PHI/ANGSPC) where PHI is the interior angle at (XEYE,YEYE).
!
!    Output, integer NVSVRT, (number of visible vertices) - 1.
!
!    Output, double precision THETA(0:NVSVRT), the polar angles of visible 
!    vertices with respect to (XEYE,YEYE) at origin and (XC(0),YC(0))
!    on positive x-axis.
!
  integer maxn
!
  double precision alpha
  double precision ang
  double precision ang1
  double precision ang2
  double precision angdif
  double precision angle
  double precision angsp2
  double precision angspc
  double precision cosang
  integer cur
  double precision diff
  double precision dx
  double precision dy
  integer i
  integer ind
  integer ivis(0:maxn)
  integer k
  integer lr
  integer lrline
  integer n
  double precision numer
  integer nvis
  integer nvsvrt
  double precision r
  double precision sinang
  double precision theta(0:maxn)
  double precision tol
  integer top
  double precision xc(0:maxn)
  double precision xeye
  double precision yc(0:maxn)
  double precision yeye
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Shift input vertices right, and possibly remove first and last
!  vertices due to collinearity with eyepoint.
!
  angsp2 = 2.0D+00 * angspc
  cur = maxn + 1
  n = maxn

  do i = nvis, 0, -1
    cur = cur - 1
    xc(cur) = xc(i)
    yc(cur) = yc(i)
    ivis(cur) = ivis(i)
  end do

  lr = lrline ( xc(cur+1), yc(cur+1), xeye, yeye, xc(cur), yc(cur), 0.0D+00 )

  if ( lr >= 0 ) then
    cur = cur + 1
    xc(0) = xc(cur)
    yc(0) = yc(cur)
    ivis(0) = ivis(cur)
  end if

  lr = lrline ( xc(n-1), yc(n-1), xeye, yeye, xc(n), yc(n), 0.0D+00 )

  if ( lr <= 0 ) then
    n = n - 1
  end if

  alpha = atan2 ( yc(0)-yeye, xc(0)-xeye )
  ang2 = 0.0D+00
  theta(0) = 0.0D+00
  top = 0
  cur = cur + 1
!
!  Process edge from vertices of indices CUR-1, CUR.
!
  do

    ang1 = ang2
    ang2 = angle ( xc(cur), yc(cur), xeye, yeye, xc(0), yc(0) )
    angdif = ang2 - ang1

    if ( angdif <= tol ) then
 
      diff = ( ( xc(cur) - xeye )**2 + ( yc(cur) - yeye)**2 ) - &
             ( ( xc(cur-1) - xeye )**2 + ( yc(cur-1) - yeye )**2 )

      if ( diff < 0.0D+00 ) then
        xc(top) = xc(cur)
        yc(top) = yc(cur)
        ivis(top) = ivis(cur)
        theta(top) = ang2
      end if

    else

      if ( angdif >= angsp2 ) then

        k = int ( angdif / angspc )
        ind = -abs ( ivis(cur))
        angdif = angdif / dble(k)
        dx = xc(cur) - xc(cur-1)
        dy = yc(cur) - yc(cur-1)
        numer = ( xc(cur) - xeye ) * dy - ( yc(cur) - yeye ) * dx

        do i = 1, k-1
          top = top + 1
          theta(top) = ang1 + dble ( i ) * angdif
          ang = theta(top) + alpha
          cosang = cos(ang)
          sinang = sin(ang)
          r = numer / ( dy * cosang - dx * sinang )
          xc(top) = r * cosang + xeye
          yc(top) = r * sinang + yeye
          ivis(top) = ind
        end do
 
      end if

      top = top + 1
      xc(top) = xc(cur)
      yc(top) = yc(cur)
      ivis(top) = ivis(cur)
      theta(top) = ang2

    end if

    cur = cur + 1

    if ( cur > n ) then
      exit
    end if

  end do

  nvsvrt = top

  return
end
subroutine vornbr ( xeye, yeye, nvrt, xc, yc, nvor, ivor, xvor, yvor, ierror )
!
!*******************************************************************************
!
!! VORNBR determines the Voronoi neighbors of an eyepoint.
!
!
!  Purpose: 
!
!    Determine the Voronoi neighbors of (XEYE,YEYE) from a
!    list of vertices which are in increasing "polar angle" order.
!    The Voronoi neighbors are a sublist of this list.  The
!    Voronoi polygon is restricted to the sector formed from the
!    the edges joining (XEYE,YEYE) to the first and last vertices
!    of this list.  Each Voronoi neighbor corresponds to an edge
!    of the Voronoi polygon.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XEYE, YEYE, the coordinates of eyepoint.
!
!    Input, integer NVRT, (number of vertices in list) minus 1.
!
!    Input, double precision XC(0:NVRT), YC(0:NVRT), the vertex 
!    coordinates from which Voronoi neighbors are determined; (XC(0),YC(0)),...,
!    (XC(NVRT),YC(NVRT)) are in increasing angular
!    displacement order with respect to (XEYE,YEYE).
!
!    Output, integer NVOR, (number of Voronoi neighbors) minus 1 [<= NVRT].
!
!    Output, integer IVOR(0:NVOR), the indices of Voronoi neighbors in XC, YC
!    arrays; 0 <= IVOR(0) < ... < IVOR(NVOR) <= NVRT.
!
!    Workspace, double precision XVOR(0:NVRT), YVOR(0:NVRT), arrays for
!    storing the vertex coordinates of the Voronoi polygon.
!
!    Output, integer IERROR, set to 212 if an error occurred.
!
  integer nvrt
!
  double precision a11
  double precision a12
  double precision a21
  double precision a22
  double precision b1
  double precision b2
  double precision det
  integer ierror
  integer im
  integer ivor(0:nvrt)
  integer k
  integer lr
  integer lrline
  integer m
  integer nvor
  double precision tol
  double precision tolabs
  double precision xc(0:nvrt)
  double precision xeye
  double precision xi
  double precision xvor(0:nvrt)
  double precision yc(0:nvrt)
  double precision yeye
  double precision yi
  double precision yvor(0:nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  k = 1
  m = 0
  ivor(0) = 0
  xvor(0) = ( xeye + xc(0) ) * 0.5D+00
  yvor(0) = ( yeye + yc(0) ) * 0.5D+00
!
!  Beginning of main loop
!
  do while ( k <= nvrt )
!
!  Determine the intersection of the perpendicular bisectors
!  of edges from (XEYE,YEYE) to (XC(K),YC(K)) and from
!  (XEYE,YEYE) to (XC(IM),YC(IM)).
!
     im = ivor(m)

     a11 = xc(k) - xeye
     a12 = yc(k) - yeye
     a21 = xc(im) - xeye
     a22 = yc(im) - yeye

     tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
     det = a11 * a22 - a21 * a12

     if ( abs ( det ) <= tolabs ) then
       ierror = 212
       return
     end if

     b1 = ( a11**2 + a12**2 ) * 0.5D+00
     b2 = ( a21**2 + a22**2 ) * 0.5D+00

     xi = ( b1 * a22 - b2 * a12 ) / det
     yi = ( b2 * a11 - b1 * a21 ) / det
!
!  Determine whether (XVOR(M+1),YVOR(M+1)) is to the left of or
!  on the directed line from (XEYE,YEYE) to (XVOR(M),YVOR(M)).
!
     xvor(m+1) = xi + xeye
     yvor(m+1) = yi + yeye
     lr = lrline ( xvor(m+1), yvor(m+1), xeye, yeye, xvor(m), yvor(m), 0.0D+00 )

     if ( lr <= 0 ) then
       m = m + 1
       ivor(m) = k
       k = k + 1
     else if ( m > 0 ) then
       m = m - 1
     else
!
!  Determine the intersection of edge from (XEYE,YEYE) to
!  (XC(0),YC(0)) and the perpendicular bisector of the edge
!  from (XEYE,YEYE) to (XC(K),YC(K)).
!
      a11 = xc(k) - xeye
      a12 = yc(k) - yeye
      a21 = yc(0) - yeye
      a22 = xeye - xc(0)
      tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
      det = a11 * a22 - a21 * a12

      if ( abs ( det) <= tolabs ) then
        ierror = 212
        return
      end if

      b1 = ( a11**2 + a12**2 ) * 0.5D+00
      b2 = 0.0D+00
      xi = ( b1 * a22 - b2 * a12 ) / det
      yi = ( b2 * a11 - b1 * a21 ) / det
      xvor(m) = xi + xeye
      yvor(m) = yi + yeye
      ivor(m) = k
      k = k + 1

    end if

  end do
!
!  The following short loop determines which vertices at the end
!  of list are not Voronoi neighbors.
!
  do

    lr = lrline ( xvor(m), yvor(m), xeye, yeye, xc(nvrt), yc(nvrt), 0.0D+00 )

    if ( lr >= 0 ) then
      exit
    end if

    m = m - 1
    if ( m < 0 ) then
      exit
    end if

  end do

  nvor = m

  return
end
subroutine vpleft ( xc, yc, ivis )
!
!*******************************************************************************
!
!! VPLEFT is called by routine VISPOL for the LEFT operation (OPER = 1).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!     Input and updated parameters:
!
!        XC,YC,IVIS - see comments in routine VISPOL
!
  logical beye
  integer cur
  logical intsct
  integer ivis(0:*)
  integer j
  integer lr
  integer lr1
  integer lr2
  integer lrline
  integer nv
  integer oper
  integer top
  double precision xc(0:*)
  double precision xe
  double precision xu
  double precision xw
  double precision yc(0:*)
  double precision ye
  double precision yu
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a left turn, S(TOP) = V(CUR), TOP <= CUR,
!  S(TOP-1) = V(CUR-1) or on interior of edge V(CUR-1)-V(CUR).
!
10 continue

  if ( cur == nv ) then
    oper = 7
    return
  end if

  if ( .not. beye .and. top <= 2 ) then
    go to 20
  end if
!
!  Check if angular displacement of stack chain >= 2*PI or
!  interior angle at boundary viewpoint.
!
  call xedge ( 1, xe, ye, xc(nv), yc(nv), xc(top-1), yc(top-1), xc(top), &
    yc(top), xu, yu, intsct )

  if ( intsct ) then

    oper = 4
    xw = xc(cur)
    yw = yc(cur)
    lr = lrline ( xc(top), yc(top), xe, ye, xc(nv), yc(nv), 0.0D+00 )

    if ( lr == -1 ) then
      xc(top) = xu
      yc(top) = yu
      ivis(top) = -cur
    end if

    return

  end if
!
!  Process next edge.
!
20 continue

  lr = lrline ( xc(cur+1), yc(cur+1), xe, ye, xc(cur), yc(cur), 0.0D+00 )

  if ( lr == -1 ) then

    cur = cur + 1
    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = cur

  else

    j = cur + 1
    lr1 = lrline ( xc(j), yc(j), xc(top-1), yc(top-1), xc(cur), yc(cur), &
      0.0D+00 )

    if ( lr1 == 1 ) then

      oper = 3
      cur = j

    else

      if ( lr == 1 ) then
        lr2 = 1
        go to 40
      end if

      do

        j = j + 1
        lr2 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )

        if ( lr2 /= 0 ) then
          exit
        end if

      end do

40    continue

      if ( lr2 == -1 ) then
        top = top + 1
        xc(top) = xc(j-1)
        yc(top) = yc(j-1)
        ivis(top) = j - 1 + nv
        top = top + 1
        xc(top) = xc(j)
        yc(top) = yc(j)
        ivis(top) = j
      else
        oper = 2
      end if

      cur = j

    end if

  end if
!
!  This test avoids extra subroutine calls.
!
  if ( oper == 1 ) then
    go to 10
  end if

  return
end
subroutine vprght ( xc, yc, ivis, ierror )
!
!*******************************************************************************
!
!! VPRGHT is called by routine VISPOL for the RIGHT operation (OPER = 2).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!     Input and updated parameters:
!        XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 206.
!
  logical beye
  integer case
  integer cur
  integer ierror
  logical intsct
  integer ivis(0:*)
  integer j
  integer lr
  integer lr1
  integer lr2
  integer lrline
  integer nv
  integer oper
  integer top
  double precision xc(0:*)
  double precision xe
  double precision xu
  double precision xw
  double precision yc(0:*)
  double precision ye
  double precision yu
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn, EYE-S(TOP)-V(CUR) is a right
!  turn, EYE-S(TOP-1)-S(TOP) is a left turn, TOP < CUR, S(TOP) =
!  V(CUR-1) and S(TOP-1)-S(TOP)-V(CUR) is a left turn or S(TOP) is
!  not on edge V(CUR-1)-V(CUR) and V(CUR-1)-V(CUR) intersects
!  EYE-S(TOP).
!  Pop points from stack. If BEYE, it is not possible for
!  (XC(CUR),YC(CUR)) to be identical to any stack points.
!
  ierror = 0

10 continue

  case = 0
  j = top

20 continue

  if ( abs ( ivis(j)) <= nv ) then

    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(j-1), yc(j-1), 0.0D+00 )

    if ( lr == -1 ) then

      case = 1

    else if ( lr == 0 ) then

      if ( abs ( ivis(j-1)) <= nv ) then
        j = j - 1
        case = 2
      else if ( (xc(j-2) - xe)**2 + (yc(j-2) - ye)**2 >= &
                (xc(j-1) - xe)**2 + (yc(j-1) - ye)**2 ) then
        j = j - 2
        case = 2
      else
        case = -1
      end if

    end if

  else if ( case == -1 ) then

    if ( ( xc(j-1) - xe )**2 + ( yc(j-1) - ye )**2 >= &
         ( xc(cur) - xe )**2 + ( yc(cur) - ye )**2 ) then
      j = j - 1
      case = 2
    else
      xw = xc(cur)
      yw = yc(cur)
      case = 3
    end if

  else

    call xedge ( 0, xc(cur-1), yc(cur-1), xc(cur), yc(cur), &
      xc(j-1), yc(j-1), xc(j), yc(j), xw, yw, intsct )

    if ( intsct ) then
      case = 3
    end if

  end if

  if ( case > 0 ) then
    go to 30
  end if

  j = j - 1
  if ( j >= 1 ) then
    go to 20
  end if
!
!  Error from no more edges in stack.
!
  ierror = 206
  return
!
!  Process next edge.
!
30 continue

  if ( case == 3 ) then

    oper = 6
    top = j - 1

  else

    top = j
    xw = xc(cur-1)
    yw = yc(cur-1)

    if ( case == 1 ) then
      call xedge ( 1, xe, ye, xc(cur), yc(cur), xc(top-1), yc(top-1), &
            xc(top), yc(top), xu, yu, intsct )
      xc(top) = xu
      yc(top) = yu
      ivis(top) = -abs ( ivis(top))
    end if

    lr = lrline ( xc(cur+1), yc(cur+1), xe, ye, xc(cur), yc(cur), 0.0D+00 )

    if ( lr == 1 ) then

      cur = cur + 1

    else

      j = cur + 1
      lr1 = lrline ( xc(j), yc(j), xw, yw, xc(cur), yc(cur), 0.0D+00 )

      if ( lr1 == -1 ) then

        oper = 5
        cur = j

      else

        if ( lr == -1 ) then
          lr2 = -1
          go to 50
        end if

        do

          j = j + 1
          lr2 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )

          if ( lr2 /= 0 ) then
            exit
          end if

        end do

50      continue

        if ( lr2 == -1 ) then
          oper = 1
          top = top + 1
          xc(top) = xc(j-1)
          yc(top) = yc(j-1)
          ivis(top) = j - 1 + nv
          top = top + 1
          xc(top) = xc(j)
          yc(top) = yc(j)
          ivis(top) = j
        end if

        cur = j

      end if

    end if

  end if
!
!  This test avoids extra subroutine calls.
!
  if ( oper == 2 ) then
    go to 10
  end if

  return
end
subroutine vpscna ( xc, yc, ivis, ierror )
!
!*******************************************************************************
!
!! VPSCNA is called by routine VISPOL for the SCANA operation (OPER = 3).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input and updated parameters:
!    XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 207.
!
  logical beye
  integer case
  integer cur
  integer ierror
  logical intsct
  integer ivis(0:*)
  integer j
  integer k
  integer lr
  integer lr1
  integer lr2
  integer lr3
  integer lrline
  integer nv
  integer oper
  integer top
  double precision xc(0:*)
  double precision xe
  double precision xw
  double precision yc(0:*)
  double precision ye
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn or forward move, S(TOP) =
!  V(CUR-1) or EYE-S(TOP)-V(CUR-1) is a forward move and TOP = 0,
!  TOP < CUR; S(TOP-1)-S(TOP)-V(CUR) is a right turn if TOP >= 1
!  or EYE-S(TOP)-V(CUR) is a right turn if TOP = 0.
!  If BEYE, it is possible that (XC(TOP),YC(TOP)) is a non-simple
!  vertex but any edge incident on this vertex encountered during
!  scan must be invisible from (XE,YE).
!
  ierror = 0
  k = cur

10 continue

  if ( xc(k+1) == xc(top) .and. yc(k+1)  ==  yc(top) ) then

    k = k + 2

  else

    call xedge ( 1, xe, ye, xc(top), yc(top), xc(k), yc(k), xc(k+1), &
      yc(k+1), xw, yw, intsct )

    if ( intsct ) then

      lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(k), yc(k), 0.0D+00 )

      if ( lr == 1 ) then

        if ( (xc(top) - xe)**2 + (yc(top) - ye)**2 >= &
          (xw - xe)**2 + (yw - ye)**2 ) then

          if ( top > 0 ) then
            case = 1
            go to 20
          end if

        else

          lr1 = lrline ( xc(k), yc(k), xe, ye, xc(top), yc(top), 0.0D+00 )

          if ( lr1 == -1 ) then
            case = 2
            go to 20
          end if

        end if

      else

        lr1 = lrline ( xc(k+1), yc(k+1), xe, ye, xc(top), yc(top), 0.0D+00 )

        if ( lr1 == -1 ) then
          case = 3
          go to 20
        end if

      end if

    end if

    k = k + 1

  end if

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 207
  return
!
!  Process current edge.
!
20 continue

  if ( case == 3 ) then

    oper = 1
    cur = k + 1
    lr = lrline ( xc(k), yc(k), xe, ye, xc(top), yc(top), 0.0D+00 )
    top = top + 1

    if ( lr == 0 ) then
      xc(top) = xc(k)
      yc(top) = yc(k)
      ivis(top) = k + nv
    else
      xc(top) = xw
      yc(top) = yw
      ivis(top) = -(k + 1 + nv)
    end if

    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = cur

  else if ( case == 1 ) then

    cur = k + 1
    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(top), yc(top), 0.0D+00 )

    if ( lr == 1 ) then

      oper = 2

    else

      j = cur + 1
      lr1 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )
      lr2 = lrline ( xc(j), yc(j), xc(k), yc(k), xc(cur), yc(cur), 0.0D+00 )

      if ( lr1 <= 0 .and. lr2 == -1 ) then

        oper = 5
        xw = xc(k)
        yw = yc(k)
        cur = j

      else

        if ( lr1 /= 0 ) then
          lr3 = lr1
          go to 40
        end if

        do

          j = j + 1
          lr3 = lrline ( xc(j), yc(j), xe, ye, xc(cur), yc(cur), 0.0D+00 )
 
          if ( lr3 /= 0 ) then
            exit
          end if

        end do

40      continue

        if ( lr3 == 1 ) then
          oper = 2
        else
          oper = 1
          top = top + 1
          xc(top) = xc(j-1)
          yc(top) = yc(j-1)
          ivis(top) = j - 1 + nv
          top = top + 1
          xc(top) = xc(j)
          yc(top) = yc(j)
          ivis(top) = j
        end if

        cur = j

      end if

    end if

  else

    oper = 6
    cur = k + 1
    lr = lrline ( xc(cur), yc(cur), xe, ye, xc(top), yc(top), 0.0D+00 )

    if ( lr == 0 ) then
      xw = xc(cur)
      yw = yc(cur)
    end if

  end if

  return
end
subroutine vpscnb ( xc, yc, ivis, ierror )
!
!*******************************************************************************
!
!! VPSCNB is called by routine VISPOL for the SCANB operation (OPER = 4).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!     Input and updated parameters:
!        XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 208.
!
  logical beye
  integer cur
  integer ierror
  logical intsct
  integer ivis(0:*)
  integer k
  integer lr
  integer lr1
  integer lrline
  integer nv
  integer oper
  double precision tol
  double precision tolabs
  integer top
  double precision xc(0:*)
  double precision xe
  double precision xu
  double precision xw
  double precision yc(0:*)
  double precision ye
  double precision yu
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  EYE-V(CUR-1)-V(CUR) is a left turn, S(TOP) = V(CUR) or S(TOP) is
!  on interior of edge V(CUR-1)-V(CUR), TOP <= CUR, S(TOP) has
!  angular displacement of 2*PI or interior angle at boundary eye.
!  (XW,YW) is the input version of (XC(CUR),YC(CUR)).
!  If BEYE, it is possible that (XC(TOP),YC(TOP)) is a non-simple
!  point but any edge containing this point encountered during scan
!  must be invisible from (XE,YE), except for 1 case where K = CUR.
!
  tolabs = tol * ( ( xc(nv) - xc(top) )**2 + ( yc(nv) - yc(top))**2 )
  k = cur

  if ( ivis(top) < 0 .or. k + 1 == nv ) then
    go to 10
  end if

  lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(top), yc(top), 0.0D+00 )

  lr1 = lrline ( xc(k+1), yc(k+1), xc(top-1), yc(top-1), xc(top), yc(top), &
    0.0D+00 )

  if ( lr == 1 .and. lr1  ==  -1 ) then
    oper = 2
    cur = k + 1
    return
  else
    k = k + 1
  end if

10 continue

  if ( k + 1 == nv ) then

    oper = 7
    cur = nv
    top = top + 1
    xc(top) = xc(nv)
    yc(top) = yc(nv)
    ivis(top) = nv
    return

  else

      if ( k == cur ) then
         call xedge ( 0, xc(nv), yc(nv), xc(top), yc(top), xw, yw, &
           xc(k+1), yc(k+1), xu, yu, intsct )
      else
         call xedge ( 0, xc(nv), yc(nv), xc(top), yc(top), xc(k), yc(k), &
           xc(k+1), yc(k+1), xu, yu, intsct )
      end if

      if ( intsct ) then
         if ( ( xc(top) - xu )**2 + ( yc(top) - yu )**2 <= tolabs ) then
           go to 20
         end if
         lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(nv), yc(nv), 0.0D+00 )
         if ( lr == 1 ) then
           oper = 2
           cur = k + 1
           return
         end if
      end if

20    continue

      k = k + 1

   end if

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 208

  return
end
subroutine vpscnc ( xc, yc, ivis, ierror )
!
!*******************************************************************************
!
!! VPSCNC is called by routine VISPOL for the SCANC operation (OPER = 5).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input/output, XC, YC, IVIS - see comments in routine VISPOL
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 209.
!
  logical beye
  integer cur
  integer ierror
  logical intsct
  integer ivis(0:*)
  integer j
  integer k
  integer lr
  integer lr1
  integer lr2
  integer lrline
  integer nv
  integer oper
  integer top
  double precision xc(0:*)
  double precision xe
  double precision xp
  double precision xu
  double precision xw
  double precision yc(0:*)
  double precision ye
  double precision yp
  double precision yu
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a left turn or forward move, EYE-V(CUR-2)-
!  V(CUR-1) is a right turn, V(CUR-2)-V(CUR-1)-V(CUR) is a left turn,
!  TOP < CUR-1, W = V(CUR-2), S(TOP) is not on V(CUR-1)-V(CUR), EYE-
!  S(TOP)-V(CUR-1) is a backward move, EYE-S(TOP-1)-S(TOP) is a left
!  turn. If BEYE, it is possible that V(CUR-1) is a non-simple point,
!  but intersection at (XC(TOP),YC(TOP)) cannot occur.
!
  ierror = 0
  xp = xc(cur-1)
  yp = yc(cur-1)
  k = cur

10 continue

  if ( xc(k+1) == xp .and. yc(k+1) == yp ) then

    go to 40

  else if ( xc(k) == xp .and. yc(k) == yp ) then

      j = k + 1
      lr = lrline ( xc(j), yc(j), xe, ye, xp, yp, 0.0D+00 )
      lr1 = lrline ( xc(j), yc(j), xw, yw, xp, yp, 0.0D+00 )

      if ( lr <= 0 .and. lr1 == -1 ) then
        go to 40
      end if

      if ( lr /= 0 ) then

        lr2 = lr

      else

        do

          j = j + 1
          lr2 = lrline ( xc(j), yc(j), xe, ye, xp, yp, 0.0D+00 )
 
          if ( lr2 /= 0 ) then
            exit
          end if

        end do

      end if

      if ( lr2 == 1 ) then
        oper = 2
      else
        oper = 1
        top = top + 1
        xc(top) = xc(j-1)
        yc(top) = yc(j-1)
        ivis(top) = j - 1 + nv
        top = top + 1
        xc(top) = xc(j)
        yc(top) = yc(j)
        ivis(top) = j
      end if

      cur = j
      return

   else

      call xedge ( 0, xp, yp, xc(top), yc(top), xc(k), yc(k), xc(k+1), &
        yc(k+1), xu, yu, intsct )

      if ( intsct ) then
        lr = lrline ( xc(k+1), yc(k+1), xe, ye, xp, yp, 0.0D+00 )
        if ( lr == 1 ) then
          oper = 2
          cur = k + 1
          return
        end if
      end if

   end if

40    continue

   k = k + 1

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 209

  return
end
subroutine vpscnd ( xc, yc, ivis, ierror )
!
!*******************************************************************************
!
!! VPSCND is called by routine VISPOL for the SCAND operation (OPER = 6).
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input/output, XC,YC,IVIS - see comments in routine VISPOL
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 210
!
  logical beye
  integer cur
  integer ierror
  logical intsct
  integer ivis(0:*)
  integer k
  integer lr
  integer lr1
  integer lr2
  integer lrline
  integer nv
  integer oper
  integer top
  double precision xc(0:*)
  double precision xe
  double precision xp
  double precision xu
  double precision xw
  double precision yc(0:*)
  double precision ye
  double precision yp
  double precision yu
  double precision yw
!
  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn, S(TOP) is a V vertex not on
!  V(CUR-1)-V(CUR), TOP < CUR, W is intersection of V(CUR-1)-V(CUR)
!  and ray EYE-S(TOP), EYE-S(TOP)-W is a forward move, and
!  EYE-S(TOP-1)-S(TOP) is a left turn if TOP >= 1.
!  If BEYE, it is possible that (XW,YW) is a non-simple point,
!  but intersection at (XC(TOP),YC(TOP)) cannot occur.
!
  ierror = 0
  xp = xc(cur-1)
  yp = yc(cur-1)
  k = cur

10 continue

  call xedge ( 0, xw, yw, xc(top), yc(top), xc(k), yc(k), xc(k+1), yc(k+1), &
    xu, yu, intsct )

  if ( intsct ) then

    lr = lrline ( xc(k+1), yc(k+1), xe, ye, xc(k), yc(k), 0.0D+00 )
    lr1 = lrline ( xc(k+1), yc(k+1), xe, ye, xc(top), yc(top), 0.0D+00 )

    if ( lr == -1 .and. lr1  ==  -1 ) then

      if ( xc(k) /= xw .or. yc(k) /= yw ) then
        go to 20
      end if

      lr2 = lrline ( xc(k+1), yc(k+1), xp, yp, xw, yw, 0.0D+00 )

      if ( lr2 == -1 ) then
        go to 30
      end if

20    continue

         oper = 1
         cur = k + 1
         lr2 = lrline ( xc(k), yc(k), xe, ye, xc(top), yc(top), 0.0D+00 )
         top = top + 1
         if ( lr2 == 0 ) then
          xc(top) = xc(k)
          yc(top) = yc(k)
          ivis(top) = k + nv
         else
          xc(top) = xu
          yc(top) = yu
          ivis(top) = -(k + 1 + nv)
         end if
         top = top + 1
         xc(top) = xc(cur)
         yc(top) = yc(cur)
         ivis(top) = cur
         return
      end if

  end if

30 continue

  k = k + 1

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierror = 210

  return
end
subroutine walkt2 ( x, y, ntri, vcl, til, tnbr, itri, iedg, ierror )
!
!*******************************************************************************
!
!! WALKT2 searches for a triangle containing a point.
!
!
!  Purpose: 
!
!    Walk through neighboring triangles of a 2-D Delaunay
!    triangulation until a triangle is found containing point (X,Y)
!    or (X,Y) is found to be outside the convex hull.  Search is
!    guaranteed to terminate for a Delaunay triangulation, else a
!    cycle may occur.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a 2-D point.
!
!    Input, integer NTRI, the number of triangles in the triangulation; used 
!    to detect cycle.
!
!    Input, double precision VCL(2,1:*), the coordinates of 2-D vertices.
!
!    Input, integer TIL(3,NTRI), the triangle incidence list.
!
!    Input, integer TNBR(3,NTRI), the triangle neighbor list.
!
!    Input/output, integer ITRI.  On input, the index of triangle to begin 
!    search at.  On output, the index of triangle that search ends at.
!
!    Output, integer IEDG, indicates the position of the point (X,Y) in
!    triangle ITRI.  A small tolerance is allowed in positions:
!    0, the interior of the triangle; 
!    1, interior of edge 1;
!    2, interior of edge 2;
!    3, interior or edge 3;
!    4, vertex 1;
!    5, vertex 2;
!    6, vertex 3;
!    -1, outside convex hull, past edge 1;
!    -2, outside convex hull, past edge 2;
!    -3, outside convex hull, past edge 3.
!
!    Output, integer IERROR, error flag.  On abnormal return,
!    IERROR is set to 226.
!
  integer ntri
!
  integer a
  double precision alpha
  integer b
  double precision beta
  integer c
  integer cnt
  double precision det
  double precision dx
  double precision dxa
  double precision dxb
  double precision dy
  double precision dya
  double precision dyb
  double precision gamma
  integer i
  integer iedg
  integer ierror
  integer itri
  integer til(3,ntri)
  integer tnbr(3,ntri)
  double precision tol
  double precision vcl(2,*)
  double precision x
  double precision y
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )

  cnt = 0
  iedg = 0
  ierror = 0

  do

    cnt = cnt + 1

    if ( cnt > ntri ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WALKT2 - Fatal error!'
      write ( *, '(a)' ) '  All triangles have been searched.'
      ierror = 226
      return
    end if
!
!  Get the vertices of triangle ITRI.
!
    a = til(1,itri)
    b = til(2,itri)
    c = til(3,itri)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point (X,Y).
!
    dxa = vcl(1,a) - vcl(1,c)
    dya = vcl(2,a) - vcl(2,c)

    dxb = vcl(1,b) - vcl(1,c)
    dyb = vcl(2,b) - vcl(2,c)

    dx = x - vcl(1,c)
    dy = y - vcl(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point (X,Y) with respect
!  to this triangle.
!
    alpha = ( dx * dyb - dy * dxb ) / det
    beta = ( dxa * dy - dya * dx ) / det
    gamma = 1.0D+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle.
!
    if ( alpha > tol .and. beta > tol .and. gamma > tol ) then
      exit
    end if
!
!  If any barycentric coordinate is (strongly) negative with respect to
!  a side, and if that side is on the convex hull, the point is outside
!  the triangles, and we are done.
!
    if ( alpha < -tol ) then
      i = tnbr(2,itri)
      if ( i <= 0 ) then
        iedg = -2
        exit
      end if
    else if ( beta < -tol ) then
      i = tnbr(3,itri)
      if ( i <= 0 ) then
        iedg = -3
        exit
      end if
    else if ( gamma < -tol ) then
      i = tnbr(1,itri)
      if ( i <= 0 ) then
        iedg = -1
        exit
      end if
!
!  At least one barycentric coordinate is between -TOL and TOL,
!  and no barycentric coordinate is less than -TOL.  We are going
!  to assign the position to an edge or vertex.
!
    else if ( alpha <= tol ) then
      if ( beta <= tol ) then
        iedg = 6
      else if ( gamma <= tol ) then
        iedg = 5
      else
        iedg = 2
      end if
      exit
    else if ( beta <= tol ) then
      if ( gamma <= tol ) then
        iedg = 4
      else
        iedg = 3
      end if
      exit
    else
      iedg = 1
      exit
    end if
!
!  If we fell through, then at least one barycentric coordinate was negative
!  for a side of the current triangle, and that side has a neighboring
!  triangle I.  Let's go there.
!
    itri = i

  end do

  return
end
subroutine width2 ( nvrt, xc, yc, i1, i2, widsq, ierror )
!
!*******************************************************************************
!
!! WIDTH2 finds the minimum breadth of a convex polygon.
!
!
!  Discussion:
!
!    WIDTH2 finds the width (minimum breadth) of a convex polygon with
!    vertices given in counter-clockwise order and with all interior 
!    angles < PI.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NVRT, the number of vertices.
!
!    Input, double precision XC(1:NVRT), YC(1:NVRT), the vertex coordinates, in
!    counter-clockwise order.
!
!    Output, integer I1, I2, indices in XC, YC such that the width is 
!    the distance from vertex (XC(I1),YC(I1)) to the line joining 
!    (XC(I2),YC(I2)) and (XC(I2+1),YC(I2+1)), where index NVRT+1 
!    is same as 1.
!
!    Output, double precision WIDSQ, the square of the width of the polygon.
!
!    Output, integer IERROR, the error flag.
!    0, no error was detected.
!    201, an error was detected.
!
  integer nvrt
!
  integer a
  double precision area1
  double precision area2
  double precision areatr
  integer b
  integer c
  double precision c1mtol
  double precision c1ptol
  double precision dist
  double precision dx
  double precision dy
  integer i_wrap
  integer i1
  integer i2
  integer ierror
  integer j
  integer jp1
  integer k
  integer kp1
  integer m
  double precision tol
  double precision widsq
  double precision xc(nvrt)
  double precision yc(nvrt)
!
  ierror = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Find the first vertex which is farthest from the edge connecting
!  vertices NVRT and 1.
!
  c1mtol = 1.0D+00 - tol
  c1ptol = 1.0D+00 + tol
  j = nvrt
  jp1 = 1
  k = 2
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

    if ( area2 <= area1 * c1ptol ) then
      exit
    end if

    area1 = area2
    k = k + 1

  end do

  m = k
  widsq = 0.0D+00
!
!  Find width = minimum distance of antipodal edge-vertex pairs.
!
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

  do

    kp1 = i_wrap ( k+1, 1, nvrt )

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(kp1), yc(kp1) )

    if ( area2 > area1 * c1ptol ) then

      a = j
      b = k
      k = k + 1
      c = k
      if ( c > nvrt ) then
        c = 1
      end if
      area1 = area2

    else if ( area2 < area1 * c1mtol ) then

      a = k
      b = j
      c = jp1
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

    else

      a = k
      b = j
      c = jp1
      k = k + 1
      j = jp1
      jp1 = j + 1
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

    end if

    if ( j > m .or. k > nvrt ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WIDTH2 - Fatal error!'
      write ( *, '(a)' ) '  J > M or K > NVRT'
      write ( *, '(a,i6)' ) '  J = ', j
      write ( *, '(a,i6)' ) '  M = ', m
      write ( *, '(a,i6)' ) '  K = ', k
      write ( *, '(a,i6)' ) '  NVRT = ', nvrt
      ierror = 201
      return
    end if

    dx = xc(c) - xc(b)
    dy = yc(c) - yc(b)
    dist = ( ( yc(a) - yc(b) ) * dx - ( xc(a) - xc(b) ) * dy )**2 &
      / ( dx**2 + dy**2 )

    if ( dist < widsq .or. widsq <= 0.0D+00 ) then
      widsq = dist
      i1 = a
      i2 = b
    end if

    if ( j == m .and. k == nvrt ) then
      exit
    end if

  end do

  return
end
subroutine xedge ( mode, xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, xu, yu, intsct )
!
!*******************************************************************************
!
!! XEDGE determines if an edge intersects another edge or ray.
!
!
!  Discussion:
!
!    An edge is a finite line segment.  A ray is a semi-infinite line
!    segment.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer MODE, is 0 for two edges, 1 (or nonzero) for a ray 
!    and an edge.
!
!    Input, double precision XV1, YV1, XV2, YV2, XW1, YW1, XW2, YW2, the
!    vertex coordinates;  an edge (ray) is from (XV1,YV1) to (thru) (XV2,YV2);
!    an edge joins vertices (XW1,YW1) and (XW2,YW2).
!
!    Output, double precision XU, YU, the coordinates of the point of 
!    intersection iff INTSCT is .TRUE.
!
!    Output, logical INTSCT, .TRUE. if the edges/ray are nondegenerate, not
!    parallel, and intersect, .FALSE. otherwise.
!
  double precision denom
  double precision dxv
  double precision dxw
  double precision dyv
  double precision dyw
  logical intsct
  integer mode
  double precision t
  double precision tol
  double precision tolabs
  double precision xu
  double precision xv1
  double precision xv2
  double precision xw1
  double precision xw2
  double precision yu
  double precision yv1
  double precision yv2
  double precision yw1
  double precision yw2
!
  tol = 100.0D+00 * epsilon ( tol )

  intsct = .false.
  dxv = xv2 - xv1
  dyv = yv2 - yv1
  dxw = xw2 - xw1
  dyw = yw2 - yw1
  tolabs = tol * max ( abs ( dxv ), abs ( dyv ), abs ( dxw ), abs ( dyw ) )
  denom = dyv * dxw - dxv * dyw

  if ( abs ( denom ) <= tolabs ) then
    return
  end if

  t = ( dyv * ( xv1 - xw1 ) - dxv * ( yv1 - yw1 ) ) / denom

  if ( t < -tol .or. t > 1.0D+00 + tol ) then
    return
  end if

  xu = xw1 + t * dxw
  yu = yw1 + t * dyw

  if ( abs ( dxv ) >= abs ( dyv ) ) then
    t = ( xu - xv1 ) / dxv
  else
    t = ( yu - yv1 ) / dyv
  end if

  if ( mode == 0 ) then
    if ( t >= -tol .and. t <= 1.0D+00 + tol ) then
      intsct = .true.
    end if
  else
    if ( t >= -tol ) then
      intsct = .true.
    end if
  end if

  return
end
subroutine xline ( xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, dv, dw, xu, &
  yu, parall )
!
!*******************************************************************************
!
!! XLINE finds the intersection of lines parallel to two other lines.
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XV1, YV1, XV2, YV2, XW1, YW1, XW2, YW2, the 
!    vertex coordinates; the first line is parallel to and at signed distance 
!    DV to the left of directed line from (XV1,YV1) to (XV2,YV2);
!    second line is parallel to and at signed distance DW to
!    left of directed line from (XW1,YW1) to (XW2,YW2)
!
!    Input, double precision DV, DW, the signed distances (positive for left).
!
!    Output, double precision XU, YU, the coordinates of the point of
!    intersection, if PARALL is .FALSE.
!
!    Output, logical PARALL, is .TRUE. if the lines are parallel, or two 
!    points for a line are identical, .FALSE. otherwise.
!
  double precision a11
  double precision a12
  double precision a21
  double precision a22
  double precision b1
  double precision b2
  double precision det
  double precision dv
  double precision dw
  logical parall
  double precision tol
  double precision tolabs
  double precision xu
  double precision xv1
  double precision xv2
  double precision xw1
  double precision xw2
  double precision yu
  double precision yv1
  double precision yv2
  double precision yw1
  double precision yw2
!
  tol = 100.0D+00 * epsilon ( tol )

  parall = .true.

  a11 = yv2 - yv1
  a12 = xv1 - xv2
  a21 = yw2 - yw1
  a22 = xw1 - xw2
  tolabs = tol * max ( abs ( a11), abs ( a12), abs ( a21), abs ( a22) )
  det = a11 * a22 - a21 * a12

  if ( abs ( det ) <= tolabs ) then
    return
  end if

  b1 = xv1 * a11 + yv1 * a12 - dv * sqrt ( a11**2 + a12**2 )
  b2 = xw1 * a21 + yw1 * a22 - dw * sqrt ( a21**2 + a22**2 )

  xu = ( b1 * a22 - b2 * a12 ) / det
  yu = ( b2 * a11 - b1 * a21 ) / det

  parall = .false.

  return
end

