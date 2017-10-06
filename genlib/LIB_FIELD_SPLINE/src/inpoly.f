
      subroutine inpoly(nd, xd,yd,np,xp,yp,ind)
C----------------------------------------------------------------------
C     This subroutine determines whether or not an (i,j) matrix 
C     element is inside a (closed) set of (i,j)'s AND has the desired
C     cloud ID.
C
C     This is essentially the same problem as -- given a point in
C     2D, are we inside a polygon. 
C----------------------------------------------------------------------

      integer np      		! # of points in polygon
      integer nd		! # points to check 
      real xd(nd) 		! 2d-locations to check
      real yd(nd)
      real    xp(np)		! 2d-locations of polygon
      real    yp(np)		
      integer ind(nd)		! THE ANSWER : ind(i)=1 if point xd(i),yd(i) is 
                                !  in polygon 0 otherwise 
      integer in
      do j = 1,nd

            call  inpoly2(xd(j),yd(j),np,xp,yp, in) 
	       ind(j)=in

      enddo

      return 
      end

      subroutine inpoly2(xpnt,ypnt,np,xp,yp,in)
C
      parameter (pi=3.14159265358979,ttpi=2.*pi)
C
      dimension xp(np),yp(np)
      real xpnt, ypnt
      integer in
C
C----------------------------------------------------------------------
C
C  THE VALUE OF THIS FUNCTION IS NON-ZERO IF AND ONLY IF (XPNT,YPNT) 
C  IS INSIDE OR *ON* THE POLYGON DEFINED BY THE POINTS (XP(I),YP(I)), 
C  FOR I FROM 1 TO NP.
C
C  THE INPUT POLYGON DOES NOT HAVE TO BE A CLOSED POLYGON, THAT IS
C  IT DOES NOT HAVE TO BE THAT (XP(1),YP(1) = (XP(NP,YP(NP)).
C
C----------------------------------------------------------------------
C
C  DETERMINE THE NUMBER OF POINTS TO LOOK AT (DEPENDING ON WHETHER THE
C  CALLER MADE THE LAST POINT A DUPLICATE OF THE FIRST OR NOT).
C
      if (xp(np).eq.xp(1) .and. yp(np).eq.yp(1)) then
        npts = np-1
      else
        npts = np
      end if

      in = 0		! ASSUME POINT IS OUTSIDE

C --- ------------------------------------------------------------------
C --- CHECK TO SEE IF THE POINT IS ON THE POLYGON.
C --- ------------------------------------------------------------------

      do ipnt = 1,npts
         if (xpnt .eq. xp(ipnt) .and. ypnt .eq. yp(ipnt) ) then
	    in = 1
	    goto 999		! EARLY EXIT
         endif
      enddo

C --- ------------------------------------------------------------------
C --- COMPUTE THE TOTAL ANGULAR CHANGE DESCRIBED BY A RAY EMANATING 
C --- FROM THE POINT (XPNT,YPNT) AND PASSING THROUGH A POINT THAT 
C --- MOVES AROUND THE POLYGON.
C --- ------------------------------------------------------------------

      anch = 0.

      inxt = npts
      xnxt = xp(npts)
      ynxt = yp(npts)
      anxt = atan2(ynxt-ypnt, xnxt-xpnt)

      do 100 ipnt=1,npts

         ilst = inxt
         xlst = xnxt
         ylst = ynxt
         alst = anxt

         inxt = ipnt
         xnxt = xp(inxt)
         ynxt = yp(inxt)
         anxt = atan2(ynxt-ypnt, xnxt-xpnt)

         adif = anxt-alst

         if (abs(adif) .gt. pi) adif = adif - sign(ttpi,adif)

         anch = anch + adif

 100  continue

C --- ------------------------------------------------------------------
C --- IF THE POINT IS OUTSIDE THE POLYGON, THE TOTAL ANGULAR CHANGE 
C --- SHOULD BE EXACTLY ZERO, WHILE IF THE POINT IS INSIDE THE POLYGON,
C --- THE TOTAL ANGULAR CHANGE SHOULD BE EXACTLY PLUS OR MINUS TWO PI.
C --- WE JUST TEST FOR THE ABSOLUTE VALUE OF THE CHANGE BEING LESS 
C --- THAN OR EQUAL TO PI.
C --- ------------------------------------------------------------------

      if (abs(anch) .ge. pi) in = 1

 999  continue
      return
      end
