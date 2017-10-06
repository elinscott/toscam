module movies
use sorting
use StringManip
contains
     
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
      subroutine displayInteractive 
      implicit none
      integer npmax
      parameter ( npmax = 300000 )
      integer ind(npmax)
      real xp(npmax) , yp(npmax) , zp(npmax) , dp(npmax)
      integer np , ic1 , id
      real psu(11) , aexpn, F(8,3)
      integer P
      real SIGN, CONTRA, BRIGHT
      real Phi, Theta
      logical Always
!***********************************************
      np=1000
      write(*,*) 'np = ', np
      do ic1=1,np
        xp(ic1)=drand1()*10.
        yp(ic1)=drand1()*10.
        zp(ic1)=drand1()*10.
        dp(ic1)=drand1()*2.
      enddo
      write(*,*) 'sorting particles...'
      call SortParticles ( dp , np , ind )
!***********************************************
!.... Initialize X-window
      id = InitDevice ( '/xwin' )
!.... set up a view point
      Phi   = 0.0
      Theta = 120.0
      psu(1)  = 0.0    ! x-coordinate (3D) of the box center
      psu(2)  = 0.0    ! y-coordinate (3D) of the box center
      psu(3)  = 0.0    ! z-coordinate (3D) of the box center
      psu(4)  = 0.5    ! scale of the plot for zoom in/out
      psu(9)  = 0.0    ! x-coordinate (3D) of the viewing point
      psu(10) = 0.0    ! y-coordinate (3D) of the viewing point
      psu(11) = 2.5    ! z-coordinate (3D) of the viewing point
!.... set up a palette
      P      = 4
      SIGN   = 1.0
      CONTRA = 0.56
      BRIGHT = 0.5
      DO
        call SetViewAngle ( Phi, Theta, psu(5), psu(6), psu(7), psu(8) )
!....   plot particles
        write(*,*) 'plotting particles...'
        call pgbbuf ()
        call pgeras ()
        call PlotParticlesPoints ( xp , yp , zp , dp , ind , np ,  &
     &                             psu , P , SIGN , CONTRA , BRIGHT )
        call InitFrame ( F )
        call DrawFrame ( F , psu )
        call PutRedshift ( 2. )
        call pgebuf
        call Control ( Phi, Theta, psu , P , SIGN , CONTRA , BRIGHT )
      ENDDO

      END subroutine

!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     ------------------------------------------------------------
      subroutine SetViewAngle ( dPhi, dTheta, CTH, STH, CPH, SPH )
!     ------------------------------------------------------------
!     Sets sines and cosines of the viewing angle
!     input  : dPhi, dTheta - viewing angles in degrees
!     output : CPH, SPH, CTH, STH - cos(Phi),sin(Phi),cos(Theta),sin(Theta)
      implicit none
      real pi, d2r
      parameter ( pi  = 3.14159265358978245 )
      parameter ( d2r = pi / 180.0 ) 
      real dPhi, dTheta, Phi, Theta, CPH, SPH, CTH, STH
      Phi   = dPhi * d2r
      CPH   = cos(Phi)
      SPH   = sin(Phi)
      Theta = dTheta * d2r
      CTH   = cos(Theta)
      STH   = sin(Theta)
      return
      end subroutine
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     ---------------------------------
      subroutine InitGIFDevice ( name )
!     ---------------------------------
!     Open GIF file with name 'name'
      implicit none
      real ScreenSize
      character*15 name
      call pgbegin ( 0 , name , 1 , 1 )
      call pgpap ( 7.515 , 0.75 )  ! 640x480
!      call pgpap ( 3.0 , 1.0 )  ! 256x256
      call pgask ( .false. )
      call pgenv ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgslw ( 1 ) ! set minimum possible linewidth
      return
      end subroutine
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     ------------------------------------------
      integer function InitDevice ( name )
!     ------------------------------------------
!     Open an arbitrary PGPLOT device
!     input : name - the name of PGPLOT device
!             nch  - number of characters in the Device name
!                    (length of the string name)   
      implicit none
      integer nch
      character*(*) name
      real w,sd
      integer pgopen, istat
      external pgopen
      nch=LEN_TRIM(name)
      istat = pgopen ( name(1:nch) )
      if ( istat .le. 0 ) then
        write(*,*) 'error: InitDevice: cannot open device',name
        stop
      else
        InitDevice = istat
      endif
      w  = 10.0
      sd = 1.0
      call PGPAP ( w , sd )  
      call PGASK ( .false. )
      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call PGSLW ( 1 ) ! set minimum possible linewidth
      call PGSCF ( 2 )     
      return
      end function
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     --------------------------------
      subroutine DrawFrame ( F , psu )
!     --------------------------------
!     Draws a wireframe 
!    input : F(8,3)  - 3D coordinates of the 8 frame corners
!             psu(11) - plot setup
!
      implicit none
      real F(8,3), psu(11)
      real Fp(8,2) ! array of coordinates of projected corners
      integer i
      do i = 1 , 8 
        call OnScreen ( F(i,1) , F(i,2) , F(i,3) , psu ,  &
     &                  Fp(i,1) , Fp(i,2) &
     &                 )
      enddo
      call PGSCI ( 1 )  ! draw the frame with a color negative to the BG
      call Line ( Fp(1,1) , Fp(1,2) , Fp(2,1) , Fp(2,2) ) 
      call Line ( Fp(1,1) , Fp(1,2) , Fp(4,1) , Fp(4,2) ) 
      call Line ( Fp(2,1) , Fp(2,2) , Fp(3,1) , Fp(3,2) ) 
      call Line ( Fp(3,1) , Fp(3,2) , Fp(4,1) , Fp(4,2) ) 
      call Line ( Fp(5,1) , Fp(5,2) , Fp(6,1) , Fp(6,2) ) 
      call Line ( Fp(5,1) , Fp(5,2) , Fp(8,1) , Fp(8,2) ) 
      call Line ( Fp(6,1) , Fp(6,2) , Fp(7,1) , Fp(7,2) ) 
      call Line ( Fp(7,1) , Fp(7,2) , Fp(8,1) , Fp(8,2) ) 
      call Line ( Fp(1,1) , Fp(1,2) , Fp(5,1) , Fp(5,2) ) 
      call Line ( Fp(2,1) , Fp(2,2) , Fp(6,1) , Fp(6,2) ) 
      call Line ( Fp(3,1) , Fp(3,2) , Fp(7,1) , Fp(7,2) ) 
      call Line ( Fp(4,1) , Fp(4,2) , Fp(8,1) , Fp(8,2) ) 
      return
      end subroutine
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     -------------------------------------
      subroutine Line ( X1 , Y1 , X2 , Y2 ) 
!     -------------------------------------
!     Draws a line between points (X1,Y1) and (X2,Y2) 
!     with clipping disabled
      implicit none
      real X1 , X2 , Y1 , Y2
      real  Xl(2) , Yl(2)
      call PGSCLP ( 0 )    ! disable clipping at the viewport edges
      Xl(1) = X1
      Yl(1) = Y1
      Xl(2) = X2
      Yl(2) = Y2      
      call PGLine ( 2 , Xl , Yl ) 
      return
      end subroutine
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     -------------------------------------------------
      subroutine OnScreen ( X , Y , Z , psu , sX , sY )
!     -------------------------------------------------
!
!     From 3D coordinates of a point calculates a 2D screen coordinates
!
!     input  : X , Y , Z  - 3D coordinates of a point 
!                           ( in the range [-0.5+Xb,0.5+Xb] by 
!                             convention, where Xb is a coordinate 
!                             of the box center )
!              psu(1)  = Xb    ! x-coordinate (3D) of the box center
!              psu(2)  = Yb    ! y-coordinate (3D) of the box center
!              psu(3)  = Zb    ! z-coordinate (3D) of the box center
!              psu(4)  = Scale ! scale of the plot for zoom in/out
!              psu(5)  = CTH   ! cosine of the viewing angle Theta
!              psu(6)  = STH   ! sine -"-
!              psu(7)  = CPH   ! cosine of the viewing angle Phi
!              psu(8)  = SPH   ! sine -"-
!              psu(9)  = Xlook ! x-coordinate (3D) of the viewing point
!              psu(10) = Ylook ! y-coordinate (3D) of the viewing point
!              psu(11) = Zlook ! z-coordinate (3D) of the viewing point
!     output : sX , sY - 2D screen coordinates 
      implicit none
      real X, Y, Z, psu(11)
      real sX , sY
      real Xd, Yd, Zd
      Xd = X*psu(5)*psu(7) + Y*psu(5)*psu(8) - Z*psu(6) - psu(1)
      Yd = Y*psu(7) - X*psu(8) - psu(2)
      Zd = X*psu(6)*psu(7) + Y*psu(6)*psu(8) + Z*psu(5) - psu(3)
      sX = psu(4) * psu(11) * Xd / abs(psu(11)-Zd) + psu(9)
      sY = psu(4) * psu(11) * Yd / abs(psu(11)-Zd) + psu(10)
      return
      end subroutine
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     --------------------------
      subroutine InitFrame ( F )
!     --------------------------
!     Sets array F(8,3) with coordinates of the box corners      
      implicit none
      real F(8,3)
      F(1,1) = -0.5
      F(1,2) = -0.5
      F(1,3) = -0.5
      F(2,1) = -0.5
      F(2,2) =  0.5
      F(2,3) = -0.5
      F(3,1) =  0.5
      F(3,2) =  0.5
      F(3,3) = -0.5
      F(4,1) =  0.5
      F(4,2) = -0.5
      F(4,3) = -0.5
      F(5,1) = -0.5
      F(5,2) = -0.5
      F(5,3) =  0.5
      F(6,1) = -0.5
      F(6,2) =  0.5
      F(6,3) =  0.5
      F(7,1) =  0.5
      F(7,2) =  0.5
      F(7,3) =  0.5
      F(8,1) =  0.5
      F(8,2) = -0.5
      F(8,3) =  0.5
      return
      end subroutine
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     --------------------------------------------
      SUBROUTINE PALETT ( TYPE , CONTRA , BRIGHT )
!     --------------------------------------------
!     Sets a "palette" of colors in the range of defined color indices 
!    This subroutine is distributed with PGPLOT in one of the demos
      INTEGER TYPE
      REAL CONTRA, BRIGHT
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL CL(5), CR(5), CG(5), CB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA CL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA CR /0.0, 0.0, 0.0, 0.3, 1.0/
      DATA CG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA CB /0.0, 0.5, 1.0, 1.0, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
     &         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/ 
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
     &         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
     &         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
     &         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- freeze
         CALL PGCTAB(CL, CR, CG, CB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END subroutine

!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!     ------------------------------------------------
      SUBROUTINE FIDDLE ( P , SIGN , CONTRA , BRIGHT )
!     ------------------------------------------------
!     Allows to interactively change the palettes, contrast, and brightness
!     (This routine is distributed with PGPLOT library in one of the 
!     example programs
      INTEGER P, IER, PGCURS
      REAL CONTRA, BRIGHT, X, Y, SIGN
      REAL X1, Y1, X2, Y2, B1, B2, C1, C2
      CHARACTER CH, CH1
      WRITE (*,*) 'Use cursor to adjust color table:'
      WRITE (*,*) ' Keys 1,2,3,4,5 select different palettes'
      WRITE (*,*) ' Key P cycles through available palettes'
      WRITE (*,*) ' Key F adjusts contrast and brightness, with'
      WRITE (*,*) '  cursor x position setting brightness [0.0 - 1.0]'
      WRITE (*,*) '   and y position setting contrast [0.0 - 10.0]'
      WRITE (*,*) '  (Hold down F key while moving cursor to change'
      WRITE (*,*) '  contrast and brightness continuously)'
      WRITE (*,*) ' Key ! resets contrast=1.0, brightness=0.5'
      WRITE (*,*) ' Key - reverses color palette'
      WRITE (*,*) ' Key X or right mouse button exits program' 
      CALL PGQWIN(X1, X2, Y1, Y2)
      B1 = 0.0
      B2 = 1.0
      C1 = 0.0
      C2 = 10.0
      CALL PGSWIN(B1, B2, C1, C2)
 10   IER = PGCURS(X, Y, CH)
      IF (CH.EQ.CHAR(0) .OR. CH.EQ.'x' .OR. CH.EQ.'X') THEN
         CALL PGSWIN(X1, X2, Y1, Y2)
         RETURN
      ELSE IF (CH.EQ.'F' .OR. CH.EQ.'f') THEN
         BRIGHT = MAX(B1, MIN(B2,X))
         CONTRA = MAX(C1, MIN(C2,Y))
      ELSE IF (CH.EQ.'C' .OR. CH.EQ.'c') THEN
         CONTRA = 1.0
         Y = 1.0
         BRIGHT = 0.5
         X = 0.5
      ELSE IF (CH.EQ.'-') THEN
         SIGN = -SIGN
      ELSE IF (CH.EQ.'1') THEN
         P = 1
      ELSE IF (CH.EQ.'2') THEN
         P = 2
      ELSE IF (CH.EQ.'3') THEN
         P = 3
      ELSE IF (CH.EQ.'4') THEN
         P = 4
      ELSE IF (CH.EQ.'5') THEN
         P = 5
      ELSE IF (CH.EQ.'P' .OR. CH.EQ.'p') THEN
         P = 1 + MOD(P,5)
      ENDIF
      CALL PALETT(P, SIGN*CONTRA, BRIGHT)
!     
      GOTO 10
      end subroutine
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!     ---------------------------------------------------------------
      subroutine Control ( Phi, Theta, psu, P, SIGN, CONTRA, BRIGHT )
!     ---------------------------------------------------------------
      implicit none
      integer IERR
      real X, Y, X1, Y1, X2, Y2, B1, B2, C1, C2
      character CH
      real Phi, Theta, psu(11)
      integer P
      real SIGN, CONTRA, BRIGHT
      integer pgcurs
      external pgcurs
      B1 = -0.5
      B2 =  0.5
      C1 = -0.5
      C2 =  0.5
      call pgswin ( B1 , B2 , C1 , C2 )
      write(*,*) 'Press R or r to rotate...'
      write(*,*) 'Press F or f to fiddle...'
      write(*,*) 'Press + or z to zoom in...'
      write(*,*) 'Press - or o to zoom out...'
      write(*,*) 'Press X or x to quit the demo...'
 10   IERR = pgcurs ( X , Y , CH )
      if ( CH .eq. CHAR(0) .or. CH.eq.'x' .or. CH.eq.'X') then
         STOP
      else if ( CH.eq.'r' .or. CH.eq.'R') then
        Theta = Theta + 10.0
        return
      else if ( CH.eq.'f' .or. CH.eq.'F') then
        call FIDDLE ( P , SIGN , CONTRA , BRIGHT )
        return
      else if ( CH.eq.'+' .or. CH.eq.'z') then
        psu(4) = psu(4) + 0.2
        return
      else if ( CH.eq.'-' .or. CH.eq.'o') then
        psu(4) = psu(4) - 0.2
        if ( psu(4) .lt. 0.0 ) psu(4) = 0.0
        return
      endif
      goto 10
      end subroutine
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!     --------------------------------
      subroutine PutRedshift ( aexpn )
!     --------------------------------
      implicit none
      integer irs , ndigits , idigit
      real aexpn, rs , frs , X , Y
      character idigits(2) , dch 
      character fdigits(2)
      character*7 time

      rs = 1.0 / aexpn - 1.0
      irs = int(rs)
      frs = rs - irs

      idigit = int(float(irs)/10)
      idigits(1) = d2ch ( idigit )
      idigit = irs - idigit*10
      idigits(2) = d2ch ( idigit )
      irs = int(frs*100)
      idigit = int(float(irs)/10)
      fdigits(1) = d2ch ( idigit )
      idigit = irs - idigit*10
      fdigits(2) = d2ch ( idigit )

      if ( idigits(1) .eq. '0' ) then
        time = 'Z='//idigits(2)//'.'//fdigits(1)//fdigits(2)
      else
       time = 'Z='//idigits(1)//idigits(2)//'.'//fdigits(1)//fdigits(2)
      endif

      call PGSCI ( 1 )
      call PGSCH (1.4)
      call PGSLW (3)
      call PGSTBG(0)
      X = -0.4
      Y = 0.5
      call PGPTXT ( X , Y , 0.0 , 0.5 , time )
      return
      end subroutine
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!     ---------------------------------------------
      subroutine SortParticles ( dp , npart , ind )
!     ---------------------------------------------
!     Sort particles according to their local density
!     ind(np) has index of particle with largest local density
!     input  : npart     - number of particles to sort
!              dp(npart) - density (or any other attribute) to use in sorting
!     output : ind(npart)- particle indices in sorted order
      implicit none
      integer npart
      real dp(npart)
      integer ind(npart)
      call indexx ( npart , dp , ind )
      return
      end subroutine
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!     ---------------------------------------------------------------
      subroutine PlotParticlesPoints ( xp, yp, zp, dp, ind, npart, &
     &                                 psu, P, SIGN, CONTRA, BRIGHT )
!     ---------------------------------------------------------------
!     Plots input particles coloring them according to the input densities
      implicit none
      integer npart
      integer ind(npart)
      real xp(npart), yp(npart), zp(npart), dp(npart)
      real psu(11)
      integer P
      real CONTRA , BRIGHT , SIGN
      integer CI1 , CI2 , ic1 , ip , icolor
      real dmin , dmax , dlmin , dlmax, xpr, ypr

      dmin = 1.e7
      dmax = -dmin
      do ic1 = 1 , npart
        if ( dp(ic1) .lt. dmin ) dmin = dp(ic1)
        if ( dp(ic1) .gt. dmax ) dmax = dp(ic1)
      enddo
      dlmax = log(dmax)
      if ( dmin .gt. 0.0 ) then
        dlmin = log(dmin)
      else
        dlmin = -3
      endif

      CI1 = 16
      CI2 = 86
      call PGSCIR ( CI1 , CI2 )
      call Palett ( P , SIGN*CONTRA , BRIGHT )
      call PGSITF ( 0 )
      call PGSCLP ( 0 )  ! no clipping

      do ic1 = 1 , npart
        ip = ind(ic1)
        call OnScreen ( xp(ip) , yp(ip) , zp(ip) , psu , xpr , ypr )
        if ( dp(ip) .gt. 0.0 ) then
          if ( dp(ip).gt.dmax ) then
            icolor = 1
          else
            icolor = int((CI2-CI1)*(log(dp(ip))-dlmin)/(dlmax-dlmin)) &
     &                + CI1
          endif
        else
          icolor = CI1
        endif
        call pgsci ( icolor )
        call pgpt1 ( xpr , ypr , 1 )
      enddo
      return
      end subroutine
!******************************************************
!******************************************************
!******************************************************
!******************************************************
!******************************************************

      end module

