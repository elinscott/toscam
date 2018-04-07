module tools_fit

! contains
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!       subroutine curv1 (n,x,y,slp1,slpn,islpsw,yp,temp,                 &
!      &                  sigma,ierr)                                     
! !                                                                       
!       integer n,islpsw,ierr 
!       real x(n),y(n),slp1,slpn,yp(n),temp(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute an interpolatory spline under tension through                 
! ! a sequence of functional values. the slopes at the two                
! ! ends of the curve may be specified or omitted.  for actual            
! ! computation of points on the curve it is necessary to call            
! ! the function curv2.                                                   
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of values to be interpolated (n.ge.2).              
! !                                                                       
! !   x is an array of the n increasing abscissae of the                  
! !   functional values.                                                  
! !                                                                       
! !   y is an array of the n ordinates of the values, (i. e.              
! !   y(k) is the functional value corresponding to x(k) ).               
! !                                                                       
! !   slp1 and slpn contain the desired values for the first              
! !   derivative of the curve at x(1) and x(n), respectively.             
! !   the user may omit values for either or both of these                
! !   parameters and signal this with islpsw.                             
! !                                                                       
! !   islpsw contains a switch indicating which slope data                
! !   should be used and which should be estimated by this                
! !   subroutine,                                                         
! !          = 0 if slp1 and slpn are to be used,                         
! !          = 1 if slp1 is to be used but not slpn,                      
! !          = 2 if slpn is to be used but not slp1,                      
! !          = 3 if both slp1 and slpn are to be estimated                
! !              internally.                                              
! !                                                                       
! !   yp is an array of length at least n.                                
! !                                                                       
! !   temp is an array of length at least n which is used for             
! !   scratch storage.                                                    
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the curviness desired. if abs(sigma) is nearly zero                 
! !   (e.g. .001) the resulting curve is approximately a                  
! !   cubic spline. if abs(sigma) is large (e.g. 50.) the                 
! !   resulting curve is nearly a polygonal line. if sigma                
! !   equals zero a cubic spline results.  a standard value               
! !   for sigma is approximately 1. in absolute value.                    
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   yp contains the values of the second derivative of the              
! !   curve at the given nodes.                                           
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if x-values are not strictly increasing.                   
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, slp1, slpn, islpsw and sigma are unaltered.                
! !                                                                       
! ! this subroutine references package modules ceez, terms,               
! ! and snhcsh.                                                           
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       nm1 = n-1 
!       np1 = n+1 
!       ierr = 0 
!       if (n .le. 1) go to 8 
!       if (x(n) .le. x(1)) go to 9 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/(x(n)-x(1)) 
! !                                                                       
! ! approximate end slopes                                                
! !                                                                       
!       if (islpsw .ge. 2) go to 1 
!       slpp1 = slp1 
!       go to 2 
!     1 delx1 = x(2)-x(1) 
!       delx2 = delx1+delx1 
!       if (n .gt. 2) delx2 = x(3)-x(1) 
!       if (delx1 .le. 0. .or. delx2 .le. delx1) go to 9 
!       call ceez (delx1,delx2,sigmap,c1,c2,c3,n) 
!       slpp1 = c1*y(1)+c2*y(2) 
!       if (n .gt. 2) slpp1 = slpp1+c3*y(3) 
!     2 if (islpsw .eq. 1 .or. islpsw .eq. 3) go to 3 
!       slppn = slpn 
!       go to 4 
!     3 delxn = x(n)-x(nm1) 
!       delxnm = delxn+delxn 
!       if (n .gt. 2) delxnm = x(n)-x(n-2) 
!       if (delxn .le. 0. .or. delxnm .le. delxn) go to 9 
!       call ceez (-delxn,-delxnm,sigmap,c1,c2,c3,n) 
!       slppn = c1*y(n)+c2*y(nm1) 
!       if (n .gt. 2) slppn = slppn+c3*y(n-2) 
! !                                                                       
! ! set up right hand side and tridiagonal system for yp and              
! ! perform forward elimination                                           
! !                                                                       
!     4 delx1 = x(2)-x(1) 
!       if (delx1 .le. 0.) go to 9 
!       dx1 = (y(2)-y(1))/delx1 
!       call terms (diag1,sdiag1,sigmap,delx1) 
!       yp(1) = (dx1-slpp1)/diag1 
!       temp(1) = sdiag1/diag1 
!       if (n .eq. 2) go to 6 
!       do 5 i = 2,nm1 
!         delx2 = x(i+1)-x(i) 
!         if (delx2 .le. 0.) go to 9 
!         dx2 = (y(i+1)-y(i))/delx2 
!         call terms (diag2,sdiag2,sigmap,delx2) 
!         diag = diag1+diag2-sdiag1*temp(i-1) 
!         yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag 
!         temp(i) = sdiag2/diag 
!         dx1 = dx2 
!         diag1 = diag2 
!     5   sdiag1 = sdiag2 
!     6 diag = diag1-sdiag1*temp(nm1) 
!       yp(n) = (slppn-dx1-sdiag1*yp(nm1))/diag 
! !                                                                       
! ! perform back substitution                                             
! !                                                                       
!       do 7 i = 2,n 
!         ibak = np1-i 
!     7   yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1) 
!       return 
! !                                                                       
! ! too few points                                                        
! !                                                                       
!     8 ierr = 1 
!       return 
! !                                                                       
! ! x-values not strictly increasing                                      
! !                                                                       
!     9 ierr = 2 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!       subroutine curvs (n,x,y,d,isw,s,eps,ys,ysp,sigma,temp,            &
!      &                  ierr)                                           
! !                                                                       
!       integer n,isw,ierr 
!       real x(n),y(n),d(n),s,eps,ys(n),ysp(n),sigma,temp(n,9) 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a smoothing spline under tension. for a given                 
! ! increasing sequence of abscissae (x(i)), i = 1,..., n and             
! ! associated ordinates (y(i)), i = 1,..., n, the function               
! ! determined minimizes the summation from i = 1 to n-1 of               
! ! the square of the second derivative of f plus sigma                   
! ! squared times the difference of the first derivative of f             
! ! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all               
! ! functions f with two continuous derivatives such that the             
! ! summation of the square of (f(x(i))-y(i))/d(i) is less                
! ! than or equal to a given constant s, where (d(i)), i = 1,             
! ! ..., n are a given set of observation weights. the                    
! ! function determined is a spline under tension with third              
! ! derivative discontinuities at (x(i)), i = 2,..., n-1. for             
! ! actual computation of points on the curve it is necessary             
! ! to call the function curv2. the determination of the curve            
! ! is performed by subroutine curvss, the subroutine curvs               
! ! only decomposes the workspace for curvss.                             
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of values to be smoothed (n.ge.2).                  
! !                                                                       
! !   x is an array of the n increasing abscissae of the                  
! !   values to be smoothed.                                              
! !                                                                       
! !   y is an array of the n ordinates of the values to be                
! !   smoothed, (i. e. y(k) is the functional value                       
! !   corresponding to x(k) ).                                            
! !                                                                       
! !   d is a parameter containing the observation weights.                
! !   this may either be an array of length n or a scalar                 
! !   (interpreted as a constant). the value of d                         
! !   corresponding to the observation (x(k),y(k)) should                 
! !   be an approximation to the standard deviation of error.             
! !                                                                       
! !   isw contains a switch indicating whether the parameter              
! !   d is to be considered a vector or a scalar,                         
! !          = 0 if d is an array of length n,                            
! !          = 1 if d is a scalar.                                        
! !                                                                       
! !   s contains the value controlling the smoothing. this                
! !   must be non-negative. for s equal to zero, the                      
! !   subroutine does interpolation, larger values lead to                
! !   smoother funtions. if parameter d contains standard                 
! !   deviation estimates, a reasonable value for s is                    
! !   float(n).                                                           
! !                                                                       
! !   eps contains a tolerance on the relative precision to               
! !   which s is to be interpreted. this must be greater than             
! !   or equal to zero and less than or equal to one. a                   
! !   reasonable value for eps is sqrt(2./float(n)).                      
! !                                                                       
! !   ys is an array of length at least n.                                
! !                                                                       
! !   ysp is an array of length at least n.                               
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the degree to which the first derivative part of the                
! !   smoothing functional is emphasized. if sigma is nearly              
! !   zero (e. g. .001) the resulting curve is approximately a            
! !   cubic spline. if sigma is large (e. g. 50.) the                     
! !   resulting curve is nearly a polygonal line. if sigma                
! !   equals zero a cubic spline results. a standard value for            
! !   sigma is approximately 1.                                           
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   temp is an array of length at least 9*n which is used               
! !   for scratch storage.                                                
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   ys contains the smoothed ordinate values.                           
! !                                                                       
! !   ysp contains the values of the second derivative of the             
! !   smoothed curve at the given nodes.                                  
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if s is negative,                                          
! !        = 3 if eps is negative or greater than one,                    
! !        = 4 if x-values are not strictly increasing,                   
! !        = 5 if a d-value is non-positive.                              
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, d, isw, s, eps, and sigma are unaltered.                   
! !                                                                       
! ! this subroutine references package modules curvss, terms,             
! ! and snhcsh.                                                           
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! decompose temp into nine arrays and call curvss                       
! !                                                                       
!       call curvss (n,x,y,d,isw,s,eps,ys,ysp,sigma,temp(1,1),            &
!      &             temp(1,2),temp(1,3),temp(1,4),temp(1,5),             &
!      &             temp(1,6),temp(1,7),temp(1,8),temp(1,9),             &
!      &             ierr)                                                
!       return 
!       END SUBROUTINE   
! 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!                                         
!       function curv2 (t,n,x,y,yp,sigma) 
! !                                                                       
!       integer n 
!       real t,x(n),y(n),yp(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function interpolates a curve at a given point                   
! ! using a spline under tension. the subroutine curv1 should             
! ! be called earlier to determine certain necessary                      
! ! parameters.                                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a real value to be mapped onto the interpo-              
! !   lating curve.                                                       
! !                                                                       
! !   n contains the number of points which were specified to             
! !   determine the curve.                                                
! !                                                                       
! !   x and y are arrays containing the abscissae and                     
! !   ordinates, respectively, of the specified points.                   
! !                                                                       
! !   yp is an array of second derivative values of the curve             
! !   at the nodes.                                                       
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, yp, and sigma should be input                 
! ! unaltered from the output of curv1.                                   
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   curv2 contains the interpolated value.                              
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this function references package modules intrvl and                   
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       im1 = intrvl(t,x,n) 
!       i = im1+1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/(x(n)-x(1)) 
! !                                                                       
! ! set up and perform interpolation                                      
! !                                                                       
!       del1 = t-x(im1) 
!       del2 = x(i)-t 
!       dels = x(i)-x(im1) 
!       sum = (y(i)*del1+y(im1)*del2)/dels 
!       if (sigmap .ne. 0.) go to 1 
!       curv2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*                 &
!      &        (del2+dels))/(6.*dels)                                    
!       return 
!     1 sigdel = sigmap*dels 
!       call snhcsh (ss,dummy,sigdel,-1) 
!       call snhcsh (s1,dummy,sigmap*del1,-1) 
!       call snhcsh (s2,dummy,sigmap*del2,-1) 
!       curv2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/            &
!      &            (sigdel*sigmap*(1.+ss))                               
!       return 
!       END function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!       function curvd (t,n,x,y,yp,sigma) 
! !                                                                       
!       integer n 
!       real t,x(n),y(n),yp(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function differentiates a curve at a given point                 
! ! using a spline under tension. the subroutine curv1 should             
! ! be called earlier to determine certain necessary                      
! ! parameters.                                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a real value at which the derivative is to be            
! !   determined.                                                         
! !                                                                       
! !   n contains the number of points which were specified to             
! !   determine the curve.                                                
! !                                                                       
! !   x and y are arrays containing the abscissae and                     
! !   ordinates, respectively, of the specified points.                   
! !                                                                       
! !   yp is an array of second derivative values of the curve             
! !   at the nodes.                                                       
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, yp, and sigma should be input                 
! ! unaltered from the output of curv1.                                   
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   curvd contains the derivative value.                                
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this function references package modules intrvl and                   
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       im1 = intrvl(t,x,n) 
!       i = im1+1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/(x(n)-x(1)) 
! !                                                                       
! ! set up and perform differentiation                                    
! !                                                                       
!       del1 = t-x(im1) 
!       del2 = x(i)-t 
!       dels = x(i)-x(im1) 
!       sum = (y(i)-y(im1))/dels 
!       if (sigmap .ne. 0.) go to 1 
!       curvd = sum+(yp(i)*(2.*del1*del1-del2*(del1+dels))-               &
!      &             yp(im1)*(2.*del2*del2-del1*(del2+dels)))             &
!      &             /(6.*dels)                                           
!       return 
!     1 sigdel = sigmap*dels 
!       call snhcsh (ss,dummy,sigdel,-1) 
!       call snhcsh (dummy,c1,sigmap*del1,1) 
!       call snhcsh (dummy,c2,sigmap*del2,1) 
!       curvd = sum+(yp(i)*(c1-ss)-yp(im1)*(c2-ss))/                      &
!      &        (sigdel*sigmap*(1.+ss))                                   
!       return 
!       END function 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!       function curvi (xl,xu,n,x,y,yp,sigma) 
! !                                                                       
!       integer n 
!       real xl,xu,x(n),y(n),yp(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function integrates a curve specified by a spline                
! ! under tension between two given limits. the subroutine                
! ! curv1 should be called earlier to determine necessary                 
! ! parameters.                                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   xl and xu contain the upper and lower limits of inte-               
! !   gration, respectively. (sl need not be less than or                 
! !   equal to xu, curvi (xl,xu,...) .eq. -curvi (xu,xl,...) ).           
! !                                                                       
! !   n contains the number of points which were specified to             
! !   determine the curve.                                                
! !                                                                       
! !   x and y are arrays containing the abscissae and                     
! !   ordinates, respectively, of the specified points.                   
! !                                                                       
! !   yp is an array from subroutine curv1 containing                     
! !   the values of the second derivatives at the nodes.                  
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, yp, and sigma should be input                 
! ! unaltered from the output of curv1.                                   
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   curvi contains the integral value.                                  
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this function references package modules intrvl and                   
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/(x(n)-x(1)) 
! !                                                                       
! ! determine actual upper and lower bounds                               
! !                                                                       
!       xxl = xl 
!       xxu = xu 
!       ssign = 1. 
!       if (xl .lt. xu) go to 1 
!       xxl = xu 
!       xxu = xl 
!       ssign = -1. 
!       if (xl .gt. xu) go to 1 
! !                                                                       
! ! return zero if xl .eq. xu                                             
! !                                                                       
!       curvi = 0. 
!       return 
! !                                                                       
! ! search for proper intervals                                           
! !                                                                       
!     1 ilm1 = intrvl (xxl,x,n) 
!       il = ilm1+1 
!       ium1 = intrvl (xxu,x,n) 
!       iu = ium1+1 
!       if (il .eq. iu) go to 8 
! !                                                                       
! ! integrate from xxl to x(il)                                           
! !                                                                       
!       sum = 0. 
!       if (xxl .eq. x(il)) go to 3 
!       del1 = xxl-x(ilm1) 
!       del2 = x(il)-xxl 
!       dels = x(il)-x(ilm1) 
!       t1 = (del1+dels)*del2/(2.*dels) 
!       t2 = del2*del2/(2.*dels) 
!       sum = t1*y(il)+t2*y(ilm1) 
!       if (sigma .eq. 0.) go to 2 
!       call snhcsh (dummy,c1,sigmap*del1,2) 
!       call snhcsh (dummy,c2,sigmap*del2,2) 
!       call snhcsh (ss,cs,sigmap*dels,3) 
!       sum = sum+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.))            &
!      &           *yp(il)+del2*del2*(c2-ss/2.)*yp(ilm1))/                &
!      &          (sigmap*sigmap*dels*(1.+ss))                            
!       go to 3 
!     2 sum = sum-t1*t1*dels*yp(il)/6.                                    &
!      &         -t2*(del1*(del2+dels)+dels*dels)*yp(ilm1)/12.            
! !                                                                       
! ! integrate over interior intervals                                     
! !                                                                       
!     3 if (iu-il .eq. 1) go to 6 
!       ilp1 = il+1 
!       do 5 i = ilp1,ium1 
!         dels = x(i)-x(i-1) 
!         sum = sum+(y(i)+y(i-1))*dels/2. 
!         if (sigma .eq. 0.) go to 4 
!         call snhcsh (ss,cs,sigmap*dels,3) 
!         sum = sum+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/                      &
!      &            (sigmap*sigmap*(1.+ss))                               
!         go to 5 
!     4   sum = sum-(yp(i)+yp(i-1))*dels*dels*dels/24. 
!     5   continue 
! !                                                                       
! ! integrate from x(iu-1) to xxu                                         
! !                                                                       
!     6 if (xxu .eq. x(ium1)) go to 10 
!       del1 = xxu-x(ium1) 
!       del2 = x(iu)-xxu 
!       dels = x(iu)-x(ium1) 
!       t1 = del1*del1/(2.*dels) 
!       t2 = (del2+dels)*del1/(2.*dels) 
!       sum = sum+t1*y(iu)+t2*y(ium1) 
!       if (sigma .eq. 0.) go to 7 
!       call snhcsh (dummy,c1,sigmap*del1,2) 
!       call snhcsh (dummy,c2,sigmap*del2,2) 
!       call snhcsh (ss,cs,sigmap*dels,3) 
!       sum = sum+(yp(iu)*del1*del1*(c1-ss/2.)+yp(ium1)*                  &
!      &          (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.)))            &
!      &         /(sigmap*sigmap*dels*(1.+ss))                            
!       go to 10 
!     7 sum = sum-t1*(del2*(del1+dels)+dels*dels)*yp(iu)/12.              &
!      &         -t2*t2*dels*yp(ium1)/6.                                  
!       go to 10 
! !                                                                       
! ! integrate from xxl to xxu                                             
! !                                                                       
!     8 delu1 = xxu-x(ium1) 
!       delu2 = x(iu)-xxu 
!       dell1 = xxl-x(ium1) 
!       dell2 = x(iu)-xxl 
!       dels = x(iu)-x(ium1) 
!       deli = xxu-xxl 
!       t1 = (delu1+dell1)*deli/(2.*dels) 
!       t2 = (delu2+dell2)*deli/(2.*dels) 
!       sum = t1*y(iu)+t2*y(ium1) 
!       if (sigma .eq. 0.) go to 9 
!       call snhcsh (dummy,cu1,sigmap*delu1,2) 
!       call snhcsh (dummy,cu2,sigmap*delu2,2) 
!       call snhcsh (dummy,cl1,sigmap*dell1,2) 
!       call snhcsh (dummy,cl2,sigmap*dell2,2) 
!       call snhcsh (ss,dummy,sigmap*dels,-1) 
!       sum = sum+(yp(iu)*(delu1*delu1*(cu1-ss/2.)                        &
!      &            -dell1*dell1*(cl1-ss/2.))                             &
!      &          +yp(ium1)*(dell2*dell2*(cl2-ss/2.)                      &
!      &            -delu2*delu2*(cu2-ss/2.)))/                           &
!      &          (sigmap*sigmap*dels*(1.+ss))                            
!       go to 10 
!     9 sum = sum-t1*(delu2*(dels+delu1)+dell2*(dels+dell1))*             &
!      &             yp(iu)/12.                                           &
!      &         -t2*(dell1*(dels+dell2)+delu1*(dels+delu2))*             &
!      &             yp(ium1)/12.                                         
! !                                                                       
! ! correct sign and return                                               
! !                                                                       
!    10 curvi = ssign*sum 
!       return 
!       END function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
!                                         
!       subroutine curvp1 (n,x,y,p,yp,temp,sigma,ierr) 
! !                                                                       
!       integer n,ierr 
!       real x(n),y(n),p,yp(n),temp(1),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a periodic interpolatory spline under tension                 
! ! through a sequence of functional values. for actual ends              
! ! of the curve may be specified or omitted.  for actual                 
! ! computation of points on the curve it is necessary to call            
! ! the function curvp2.                                                  
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of values to be interpolated (n.ge.2).              
! !                                                                       
! !   x is an array of the n increasing abscissae of the                  
! !   functional values.                                                  
! !                                                                       
! !   y is an array of the n ordinates of the values, (i. e.              
! !   y(k) is the functional value corresponding to x(k) ).               
! !                                                                       
! !   p is the period (p .gt. x(n)-x(1)).                                 
! !                                                                       
! !   yp is an array of length at least n.                                
! !                                                                       
! !   temp is an array of length at least 2*n which is used               
! !   for scratch storage.                                                
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor.  this value indicates            
! !   the curviness desired. if abs(sigma) is nearly zero                 
! !   (e.g. .001) the resulting curve is approximately a                  
! !   cubic spline. if abs(sigma) is large (e.g. 50.) the                 
! !   resulting curve is nearly a polygonal line. if sigma                
! !   equals zero a cubic spline results.  a standard value               
! !   for sigma is approximately 1. in absolute value.                    
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   yp contains the values of the second derivative of the              
! !   curve at the given nodes.                                           
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if p is less than or equal to x(n)-x(1),                   
! !        = 3 if x-values are not strictly increasing.                   
! !                                                                       
! ! and                                                                   
! !                                                                       
! !  n, x, y, and sigma are unaltered.                                    
! !                                                                       
! ! this subroutine references package modules terms and                  
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       nm1 = n-1 
!       np1 = n+1 
!       ierr = 0 
!       if (n .le. 1) go to 6 
!       if (p .le. x(n)-x(1) .or. p .le. 0.) go to 7 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/p 
! !                                                                       
! ! set up right hand side and tridiagonal system for yp and              
! ! perform forward elimination                                           
! !                                                                       
!       delx1 = p-(x(n)-x(1)) 
!       dx1 = (y(1)-y(n))/delx1 
!       call terms (diag1,sdiag1,sigmap,delx1) 
!       delx2 = x(2)-x(1) 
!       if (delx2 .le. 0.) go to 8 
!       dx2 = (y(2)-y(1))/delx2 
!       call terms (diag2,sdiag2,sigmap,delx2) 
!       diag = diag1+diag2 
!       yp(1) = (dx2-dx1)/diag 
!       temp(np1) = -sdiag1/diag 
!       temp(1) = sdiag2/diag 
!       dx1 = dx2 
!       diag1 = diag2 
!       sdiag1 = sdiag2 
!       if (n .eq. 2) go to 2 
!       do 1 i = 2,nm1 
!         npi = n+i 
!         delx2 = x(i+1)-x(i) 
!         if (delx2 .le. 0.) go to 8 
!         dx2 = (y(i+1)-y(i))/delx2 
!         call terms (diag2,sdiag2,sigmap,delx2) 
!         diag = diag1+diag2-sdiag1*temp(i-1) 
!         yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag 
!         temp(npi) = -temp(npi-1)*sdiag1/diag 
!         temp(i) = sdiag2/diag 
!         dx1 = dx2 
!         diag1 = diag2 
!     1   sdiag1 = sdiag2 
!     2 delx2 = p-(x(n)-x(1)) 
!       dx2 = (y(1)-y(n))/delx2 
!       call terms (diag2,sdiag2,sigmap,delx2) 
!       yp(n) = dx2-dx1 
!       temp(nm1) = temp(2*n-1)-temp(nm1) 
!       if (n .eq. 2) go to 4 
! !                                                                       
! ! perform first step of back substitution                               
! !                                                                       
!       do 3 i = 3,n 
!         ibak = np1-i 
!         npibak =n+ibak 
!         yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1) 
!     3   temp(ibak) =temp(npibak)-temp(ibak)*temp(ibak+1) 
!     4 yp(n) = (yp(n)-sdiag2*yp(1)-sdiag1*yp(nm1))/                      &
!      &        (diag1+diag2+sdiag2*temp(1)+sdiag1*temp(nm1))             
! !                                                                       
! ! perform second step of back substitution                              
! !                                                                       
!       ypn =   yp(n) 
!       do 5 i = 1,nm1 
!     5    yp(i) = yp(i)+temp(i)*ypn 
!       return 
! !                                                                       
! ! too few points                                                        
! !                                                                       
!     6 ierr = 1 
!       return 
! !                                                                       
! ! period too small                                                      
! !                                                                       
!     7 ierr = 2 
!       return 
! !                                                                       
! ! x-values not strictly increasing                                      
! !                                                                       
!     8 ierr = 3 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!                                         
!       subroutine curvps (n,x,y,p,d,isw,s,eps,ys,ysp,sigma,              &
!      &                   temp,ierr)                                     
! !                                                                       
!       integer n,isw,ierr 
!       real x(n),y(n),p,d(n),s,eps,ys(n),ysp(n),sigma,                   &
!      &     temp(n,11)                                                   
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a periodic smoothing spline under tension. for a              
! ! given increasing sequence of abscissae (x(i)), i = 1,...,n            
! ! and associated ordinates (y(i)), i = 1,...,n, letting p be            
! ! the period, x(n+1) = x(1)+p, and y(n+1) = y(1), the                   
! ! function determined minimizes the summation from i = 1 to             
! ! n of the square of the second derivative of f plus sigma              
! ! squared times the difference of the first derivative of f             
! ! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all               
! ! functions f with period p and two continuous derivatives              
! ! such that the summation of the square of                              
! ! (f(x(i))-y(i))/d(i) is less than or equal to a given                  
! ! constant s, where (d(i)), i = 1,...,n are a given set of              
! ! observation weights. the function determined is a periodic            
! ! spline under tension with third derivative discontinuities            
! ! at (x(i)) i = 1,...,n (and all periodic translations of               
! ! these values). for actual computation of points on the                
! ! curve it is necessary to call the function curvp2. the                
! ! determination of the curve is performed by subroutine                 
! ! curvpp, the subroutin curvps only decomposes the workspace            
! ! for curvpp.                                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of values to be smoothed (n.ge.2).                  
! !                                                                       
! !   x is an array of the n increasing abscissae of the                  
! !   values to be smoothed.                                              
! !                                                                       
! !   y is an array of the n ordinates of the values to be                
! !   smoothed, (i. e. y(k) is the functional value                       
! !   corresponding to x(k) ).                                            
! !                                                                       
! !   p is the period (p .gt. x(n)-x(1)).                                 
! !                                                                       
! !   d is a parameter containing the observation weights.                
! !   this may either be an array of length n or a scalar                 
! !   (interpreted as a constant). the value of d                         
! !   corresponding to the observation (x(k),y(k)) should                 
! !   be an approximation to the standard deviation of error.             
! !                                                                       
! !   isw contains a switch indicating whether the parameter              
! !   d is to be considered a vector or a scalar,                         
! !          = 0 if d is an array of length n,                            
! !          = 1 if d is a scalar.                                        
! !                                                                       
! !   s contains the value controlling the smoothing. this                
! !   must be non-negative. for s equal to zero, the                      
! !   subroutine does interpolation, larger values lead to                
! !   smoother funtions. if parameter d contains standard                 
! !   deviation estimates, a reasonable value for s is                    
! !   float(n).                                                           
! !                                                                       
! !   eps contains a tolerance on the relative precision to               
! !   which s is to be interpreted. this must be greater than             
! !   or equal to zero and less than or equal to one. a                   
! !   reasonable value for eps is sqrt(2./float(n)).                      
! !                                                                       
! !   ys is an array of length at least n.                                
! !                                                                       
! !   ysp is an array of length at least n.                               
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the degree to which the first derivative part of the                
! !   smoothing functional is emphasized. if sigma is nearly              
! !   zero (e. g. .001) the resulting curve is approximately a            
! !   cubic spline. if sigma is large (e. g. 50.) the                     
! !   resulting curve is nearly a polygonal line. if sigma                
! !   equals zero a cubic spline results. a standard value for            
! !   sigma is approximately 1.                                           
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   temp is an array of length at least 11*n which is used              
! !   for scratch storage.                                                
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   ys contains the smoothed ordinate values.                           
! !                                                                       
! !   ysp contains the values of the second derivative of the             
! !   smoothed curve at the given nodes.                                  
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if s is negative,                                          
! !        = 3 if eps is negative or greater than one,                    
! !        = 4 if x-values are not strictly increasing,                   
! !        = 5 if a d-value is non-positive,                              
! !        = 6 if p is less than or equal to x(n)-x(1).                   
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, p, d, isw, s, eps, and sigma are unaltered.                
! !                                                                       
! ! this subroutine references package modules curvpp, terms,             
! ! and snhcsh.                                                           
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! decompose temp into eleven arrays and call curvpp                     
! !                                                                       
!       call curvpp (n,x,y,p,d,isw,s,eps,ys,ysp,sigma,                    &
!      &             temp(1,1),temp(1,2),temp(1,3),temp(1,4),             &
!      &             temp(1,5),temp(1,6),temp(1,7),temp(1,8),             &
!      &             temp(1,9),temp(1,10),temp(1,11),ierr)                
!       return 
!       END SUBROUTINE 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!                                           
!       function curvp2 (t,n,x,y,p,yp,sigma) 
! !                                                                       
!       integer n 
!       real t,x(n),y(n),p,yp(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function interpolates a curve at a given point using             
! ! a periodic spline under tension. the subroutine curvp1                
! ! should be called earlier to determine certain necessary               
! ! parameters.                                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a real value to be mapped onto the interpo-              
! !   lating curve.                                                       
! !                                                                       
! !   n contains the number of points which were specified to             
! !   determine the curve.                                                
! !                                                                       
! !   x and y are arrays containing the abscissae and                     
! !   ordinates, respectively, of the specified points.                   
! !                                                                       
! !   p contains the period.                                              
! !                                                                       
! !   yp is an array of second derivative values of the curve             
! !   at the nodes.                                                       
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, p, yp, and sigma should be input              
! ! unaltered from the output of curvp1.                                  
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   curvp2 contains the interpolated value.                             
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this function references package modules intrvp and                   
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       im1 = intrvp (t,x,n,p,tp) 
!       i = im1+1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/p 
! !                                                                       
! ! set up and perform interpolation                                      
! !                                                                       
!       del1 = tp-x(im1) 
!       if (im1 .eq. n) go to 1 
!       del2 = x(i)-tp 
!       dels = x(i)-x(im1) 
!       go to 2 
!     1 i = 1 
!       del2 = x(1)+p-tp 
!       dels = p-(x(n)-x(1)) 
!     2 sum = (y(i)*del1+y(im1)*del2)/dels 
!       if (sigmap .ne. 0.) go to 3 
!       curvp2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*                &
!      &        (del2+dels))/(6.*dels)                                    
!       return 
!     3 sigdel = sigmap*dels 
!       call snhcsh (ss,dummy,sigdel,-1) 
!       call snhcsh (s1,dummy,sigmap*del1,-1) 
!       call snhcsh (s2,dummy,sigmap*del2,-1) 
!       curvp2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/           &
!      &            (sigdel*sigmap*(1.+ss))                               
!       return 
!       END function
! 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
!                                      
!       function curvpi (xl,xu,n,x,y,p,yp,sigma) 
! !                                                                       
!       integer n 
!       real xl,xu,x(n),y(n),p,yp(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function integrates a curve specified by a periodic              
! ! spline under tension between two given limits. the                    
! ! subroutine curvp1 should be called earlier to determine               
! ! necessary parameters.                                                 
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   xl and xu contain the upper and lower limits of inte-               
! !   gration, respectively. (sl need not be less than or                 
! !   equal to xu, curvpi (xl,xu,...) .eq. -curvpi (xu,xl,...) ).         
! !                                                                       
! !   n contains the number of points which were specified to             
! !   determine the curve.                                                
! !                                                                       
! !   x and y are arrays containing the abscissae and                     
! !   ordinates, respectively, of the specified points.                   
! !                                                                       
! !   p contains the period.                                              
! !                                                                       
! !   yp is an array from subroutine curvp1 containing                    
! !   the values of the second derivatives at the nodes.                  
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, p, yp, and sigma should be input              
! ! unaltered from the output of curvp1.                                  
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !                                                                       
! !   curvpi contains the integral value.                                 
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this function references package modules intrvp and                   
! ! snhcsh.                                                               
! !                                                                       
! !--------------------------------------------------------------         
! !                                                                       
!       integer uper 
!       logical bdy 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/p 
! !                                                                       
! ! determine actual upper and lower bounds                               
! !                                                                       
!       x1pp = x(1)+p 
!       isign = 1 
!       ilm1 = intrvp (xl,x,n,p,xxl) 
!       lper = int((xl-x(1))/p) 
!       if (xl .lt. x(1)) lper = lper-1 
!       ium1 = intrvp (xu,x,n,p,xxu) 
!       uper = int((xu-x(1))/p) 
!       if (xu .lt. x(1)) uper = uper-1 
!       ideltp = uper-lper 
!       bdy = float(ideltp)*(xxu-xxl) .lt. 0. 
!       if ((ideltp .eq. 0 .and. xxu .lt. xxl) .or.                       &
!      &    ideltp .lt. 0) isign = -1                                     
!       if (bdy) ideltp = ideltp-isign 
!       if (xxu .ge. xxl) go to 1 
!       xsave = xxl 
!       xxl = xxu 
!       xxu = xsave 
!       isave = ilm1 
!       ilm1 = ium1 
!       ium1 = isave 
!     1 il = ilm1+1 
!       if (ilm1 .eq. n) il = 1 
!       xil = x(il) 
!       if (ilm1 .eq. n) xil = x1pp 
!       iu = ium1+1 
!       if (ium1 .eq. n) iu = 1 
!       xiu = x(iu) 
!       if (ium1 .eq. n) xiu = x1pp 
!       s1 = 0. 
!       if (ilm1 .eq. 1 .or. (ideltp .eq. 0 .and.                         &
!      &    .not. bdy)) go to 4                                           
! !                                                                       
! ! integrate from x(1) to x(ilm1), store in s1                           
! !                                                                       
!       do 3 i = 2,ilm1 
!         dels = x(i)-x(i-1) 
!         s1 = s1+(y(i)+y(i-1))*dels/2. 
!         if (sigma .eq. 0.) go to 2 
!         call snhcsh (ss,cs,sigmap*dels,3) 
!         s1 = s1+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/                        &
!      &            (sigmap*sigmap*(1.+ss))                               
!         go to 3 
!     2   s1 = s1-(yp(i)+yp(i-1))*dels*dels*dels/24. 
!     3   continue 
!     4 s2 = 0. 
!       if (x(ilm1) .ge. xxl .or. (ideltp .eq. 0                          &
!      &    .and. .not. bdy)) go to 6                                     
! !                                                                       
! ! integrate from x(ilm1) to xxl, store in s2                            
! !                                                                       
!       del1 = xxl-x(ilm1) 
!       del2 = xil-xxl 
!       dels = xil-x(ilm1) 
!       t1 = del1*del1/(2.*dels) 
!       t2 = (del2+dels)*del1/(2.*dels) 
!       s2 = t1*y(il)+t2*y(ilm1) 
!       if (sigma .eq. 0.) go to 5 
!       call snhcsh (dummy,c1,sigmap*del1,2) 
!       call snhcsh (dummy,c2,sigmap*del2,2) 
!       call snhcsh (ss,cs,sigmap*dels,3) 
!       s2 = s2+(yp(il)*del1*del1*(c1-ss/2.)+yp(ilm1)*                    &
!      &          (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.)))            &
!      &         /(sigmap*sigmap*dels*(1.+ss))                            
!       go to 6 
!     5 s2 = s2-t1*(del2*(del1+dels)                                      &
!      &     +dels*dels)*yp(il)/12.                                       &
!      &     -t2*t2*dels*yp(ilm1)/6.                                      
!     6 s3 = 0. 
!       if (xxl .ge. xil .or. (ideltp .eq. 0 .and. bdy)                   &
!      &    .or. ilm1 .eq. ium1) go to 8                                  
! !                                                                       
! ! integrate from xxl to xil, store in s3                                
! !                                                                       
!       del1 = xxl-x(ilm1) 
!       del2 = xil-xxl 
!       dels = xil-x(ilm1) 
!       t1 = (del1+dels)*del2/(2.*dels) 
!       t2 = del2*del2/(2.*dels) 
!       s3 = t1*y(il)+t2*y(ilm1) 
!       if (sigma .eq. 0.) go to 7 
!       call snhcsh (dummy,c1,sigmap*del1,2) 
!       call snhcsh (dummy,c2,sigmap*del2,2) 
!       call snhcsh (ss,cs,sigmap*dels,3) 
!       s3 = s3+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.))              &
!      &           *yp(il)+del2*del2*(c2-ss/2.)*yp(ilm1))/                &
!      &          (sigmap*sigmap*dels*(1.+ss))                            
!       go to 8 
!     7 s3 = s3-t1*t1*dels*yp(il)/6.                                      &
!      &     -t2*(del1*(del2+dels)+dels*dels)*                            &
!      &     yp(ilm1)/12.                                                 
!     8 s4 = 0. 
!       if (ilm1 .ge. ium1-1 .or. (ideltp .eq. 0 .and. bdy))              &
!      &   go to 11                                                       
! !                                                                       
! ! integrate from xil to x(ium1), store in s4                            
! !                                                                       
!       ilp1 = il+1 
!       do 10 i = ilp1,ium1 
!         dels = x(i)-x(i-1) 
!         s4 = s4+(y(i)+y(i-1))*dels/2. 
!         if (sigma .eq. 0.) go to 9 
!         call snhcsh (ss,cs,sigmap*dels,3) 
!         s4 = s4+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/                        &
!      &            (sigmap*sigmap*(1.+ss))                               
!         go to 10 
!     9   s4 = s4-(yp(i)+yp(i-1))*dels*dels*dels/24. 
!    10   continue 
!    11 s5 = 0. 
!       if (x(ium1) .ge. xxu .or. (ideltp .eq. 0 .and. bdy)               &
!      &    .or. ilm1 .eq. ium1) go to 13                                 
! !                                                                       
! ! integrate from x(ium1) to xxu, store in s5                            
! !                                                                       
!       del1 = xxu-x(ium1) 
!       del2 = xiu-xxu 
!       dels = xiu-x(ium1) 
!       t1 = del1*del1/(2.*dels) 
!       t2 = (del2+dels)*del1/(2.*dels) 
!       s5 = t1*y(iu)+t2*y(ium1) 
!       if (sigma .eq. 0.) go to 12 
!       call snhcsh (dummy,c1,sigmap*del1,2) 
!       call snhcsh (dummy,c2,sigmap*del2,2) 
!       call snhcsh (ss,cs,sigmap*dels,3) 
!       s5 = s5+(yp(iu)*del1*del1*(c1-ss/2.)+yp(ium1)*                    &
!      &          (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.)))            &
!      &         /(sigmap*sigmap*dels*(1.+ss))                            
!       go to 13 
!    12 s5 = s5-t1*(del2*(del1+dels)                                      &
!      &     +dels*dels)*yp(iu)/12.                                       &
!      &     -t2*t2*dels*yp(ium1)/6.                                      
!    13 s6 = 0. 
!       if (xxu .ge. xiu .or. (ideltp .eq. 0 .and.                        &
!      &    .not. bdy)) go to 15                                          
! !                                                                       
! ! integrate from xxu to xiu, store in s6                                
! !                                                                       
!       del1 = xxu-x(ium1) 
!       del2 = xiu-xxu 
!       dels = xiu-x(ium1) 
!       t1 = (del1+dels)*del2/(2.*dels) 
!       t2 = del2*del2/(2.*dels) 
!       s6 = t1*y(iu)+t2*y(ium1) 
!       if (sigma .eq. 0.) go to 14 
!       call snhcsh (dummy,c1,sigmap*del1,2) 
!       call snhcsh (dummy,c2,sigmap*del2,2) 
!       call snhcsh (ss,cs,sigmap*dels,3) 
!       s6 = s6+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.))              &
!      &           *yp(iu)+del2*del2*(c2-ss/2.)*yp(ium1))/                &
!      &          (sigmap*sigmap*dels*(1.+ss))                            
!       go to 15 
!    14 s6 = s6-t1*t1*dels*yp(iu)/6.                                      &
!      &     -t2*(del1*(del2+dels)+dels*dels)*                            &
!      &     yp(ium1)/12.                                                 
!    15 s7 = 0. 
!       if (iu .eq. 1 .or. (ideltp .eq. 0 .and. .not. bdy))               &
!      &    go to 18                                                      
! !                                                                       
! ! integrate from xiu to x1pp, store in s7                               
! !                                                                       
!       np1 = n+1 
!       iup1 = iu+1 
!       do 17 ii = iup1,np1 
!         im1 = ii-1 
!         i = ii 
!         if (i .eq. np1) i=1 
!         dels = x(i)-x(im1) 
!         if (dels .le. 0.) dels=dels+p 
!         s7 = s7+(y(i)+y(im1))*dels/2. 
!         if (sigma .eq. 0.) go to 16 
!         call snhcsh (ss,cs,sigmap*dels,3) 
!         s7 = s7+(yp(i)+yp(im1))*dels*(cs-ss/2.)/                        &
!      &            (sigmap*sigmap*(1.+ss))                               
!         go to 17 
!    16   s7 = s7-(yp(i)+yp(im1))*dels*dels*dels/24. 
!    17   continue 
!    18 s8 = 0. 
!       if (ilm1 .lt. ium1 .or. (ideltp .eq. 0 .and. bdy))                &
!      &   go to 20                                                       
! !                                                                       
! ! integrate from xxl to xxu, store in s8                                
! !                                                                       
!       delu1 = xxu-x(ium1) 
!       delu2 = xiu-xxu 
!       dell1 = xxl-x(ium1) 
!       dell2 = xiu-xxl 
!       dels = xiu-x(ium1) 
!       deli = xxu-xxl 
!       t1 = (delu1+dell1)*deli/(2.*dels) 
!       t2 = (delu2+dell2)*deli/(2.*dels) 
!       s8 = t1*y(iu)+t2*y(ium1) 
!       if (sigma .eq. 0.) go to 19 
!       call snhcsh (dummy,cu1,sigmap*delu1,2) 
!       call snhcsh (dummy,cu2,sigmap*delu2,2) 
!       call snhcsh (dummy,cl1,sigmap*dell1,2) 
!       call snhcsh (dummy,cl2,sigmap*dell2,2) 
!       call snhcsh (ss,dummy,sigmap*dels,-1) 
!       s8 = s8+(yp(iu)*(delu1*delu1*(cu1-ss/2.)                          &
!      &            -dell1*dell1*(cl1-ss/2.))                             &
!      &          +yp(ium1)*(dell2*dell2*(cl2-ss/2.)                      &
!      &            -delu2*delu2*(cu2-ss/2.)))/                           &
!      &          (sigmap*sigmap*dels*(1.+ss))                            
!       go to 20 
!    19 s8 = s8-t1*(delu2*(dels+delu1)                                    &
!      &     +dell2*(dels+dell1))*yp(iu)/12.                              &
!      &     -t2*(dell1*(dels+dell2)                                      &
!      &     +delu1*(dels+delu2))*yp(ium1)/12.                            
!    20 so = s1+s2+s6+s7 
!       si = s3+s4+s5+s8 
!       if (bdy) go to 21 
!       curvpi = float(ideltp)*(so+si)+float(isign)*si 
!       return 
!    21 curvpi = float(ideltp)*(so+si)+float(isign)*so 
!       return 
!       END function
! 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
! 
!                                        
!       subroutine kurv1 (n,x,y,slp1,slpn,islpsw,xp,yp,temp,s,            &
!      &                  sigma,ierr)                                     
! !                                                                       
!       integer n,islpsw,ierr 
!       real x(n),y(n),slp1,slpn,xp(n),yp(n),temp(n),s(n),                &
!      &     sigma                                                        
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a spline under tension forming a curve in the                 
! ! plane and passing through a sequence of pairs (x(1),y(1)),            
! ! ...,(x(n),y(n)). for actual computation of points on the              
! ! curve it is necessary to call the subroutine kurv2.                   
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of points to be interpolated (n.ge.2).              
! !                                                                       
! !   x is an array containing the n x-coordinates of the                 
! !   points.                                                             
! !                                                                       
! !   y is an array containing the n y-coordinates of the                 
! !   points. (adjacent x-y pairs must be distinct, i. e.                 
! !   either x(i) .ne. x(i+1) or y(i) .ne. y(i+1), for                    
! !   i = 1,...,n-1.)                                                     
! !                                                                       
! !   slp1 and slpn contain the desired values for the angles             
! !   (in radians) of the slope at (x(1),y(1)) and (x(n),y(n))            
! !   respectively. the angles are measured counter-clock-                
! !   wise from the x-axis and the positive sense of the curve            
! !   is assumed to be that moving from point 1 to point n.               
! !   the user may omit values for either or both of these                
! !   parameters and signal this with islpsw.                             
! !                                                                       
! !   islpsw contains a switch indicating which slope data                
! !   should be used and which should be estimated by this                
! !   subroutine,                                                         
! !          = 0 if slp1 and slpn are to be used,                         
! !          = 1 if slp1 is to be used but not slpn,                      
! !          = 2 if slpn is to be used but not slp1,                      
! !          = 3 if both slp1 and slpn are to be estimated                
! !              internally.                                              
! !                                                                       
! !   xp and yp are arrays of length at least n.                          
! !                                                                       
! !   temp is an array of length at least n which is used                 
! !   for scratch storage.                                                
! !                                                                       
! !   s is an array of length at least n.                                 
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the curviness desired. if abs(sigma) is nearly zero                 
! !   (e.g. .001) the resulting curve is approximately a cubic            
! !   spline. if abs(sigma) is large (e. g. 50.) the resulting            
! !   curve is nearly a polygonal line. if sigma equals zero a            
! !   cubic spline results. a standard value for sigma is                 
! !   approximately 1. in absolute value.                                 
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   xp and yp contain information about the curvature of the            
! !   curve at the given nodes.                                           
! !                                                                       
! !   s contains the polygonal arclengths of the curve.                   
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if adjacent coordinate pairs coincide.                     
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, slp1, slpn, islpsw, and sigma are unaltered.               
! !                                                                       
! ! this subroutine references package modules ceez, terms,               
! ! and snhcsh.                                                           
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       nm1 = n-1 
!       np1 = n+1 
!       ierr = 0 
!       if (n .le. 1) go to 11 
! !                                                                       
! ! determine polygonal arclengths                                        
! !                                                                       
!       s(1) = 0. 
!       do 1 i = 2,n 
!         im1 = i-1 
!     1   s(i) = s(im1)+sqrt((x(i)-x(im1))**2+                            &
!      &         (y(i)-y(im1))**2)                                        
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/s(n) 
! !                                                                       
! ! approximate end slopes                                                
! !                                                                       
!       if (islpsw .ge. 2) go to 2 
!       slpp1x = cos(slp1) 
!       slpp1y = sin(slp1) 
!       go to 4 
!     2 dels1 = s(2)-s(1) 
!       dels2 = dels1+dels1 
!       if (n .gt. 2) dels2 = s(3)-s(1) 
!       if (dels1 .eq. 0. .or. dels2 .eq. 0.) go to 12 
!       call ceez (dels1,dels2,sigmap,c1,c2,c3,n) 
!       sx = c1*x(1)+c2*x(2) 
!       sy = c1*y(1)+c2*y(2) 
!       if (n .eq. 2) go to 3 
!       sx = sx+c3*x(3) 
!       sy = sy+c3*y(3) 
!     3 delt = sqrt(sx*sx+sy*sy) 
!       slpp1x = sx/delt 
!       slpp1y = sy/delt 
!     4 if (islpsw .eq. 1 .or. islpsw .eq. 3) go to 5 
!       slppnx = cos(slpn) 
!       slppny = sin(slpn) 
!       go to 7 
!     5 delsn = s(n)-s(nm1) 
!       delsnm = delsn+delsn 
!       if (n .gt. 2) delsnm = s(n)-s(n-2) 
!       if (delsn .eq. 0. .or. delsnm .eq. 0.) go to 12 
!       call ceez (-delsn,-delsnm,sigmap,c1,c2,c3,n) 
!       sx = c1*x(n)+c2*x(nm1) 
!       sy = c1*y(n)+c2*y(nm1) 
!       if (n .eq. 2) go to 6 
!       sx = sx+c3*x(n-2) 
!       sy = sy+c3*y(n-2) 
!     6 delt = sqrt(sx*sx+sy*sy) 
!       slppnx = sx/delt 
!       slppny = sy/delt 
! !                                                                       
! ! set up right hand sides and tridiagonal system for xp and             
! ! yp and perform forward elimination                                    
! !                                                                       
!     7 dx1 = (x(2)-x(1))/s(2) 
!       dy1 = (y(2)-y(1))/s(2) 
!       call terms (diag1,sdiag1,sigmap,s(2)) 
!       xp(1) = (dx1-slpp1x)/diag1 
!       yp(1) = (dy1-slpp1y)/diag1 
!       temp(1) = sdiag1/diag1 
!       if (n .eq. 2) go to 9 
!       do 8 i = 2,nm1 
!         dels2 = s(i+1)-s(i) 
!         if (dels2 .eq. 0.) go to 12 
!         dx2 = (x(i+1)-x(i))/dels2 
!         dy2 = (y(i+1)-y(i))/dels2 
!         call terms (diag2,sdiag2,sigmap,dels2) 
!         diag = diag1+diag2-sdiag1*temp(i-1) 
!         diagin = 1./diag 
!         xp(i) = (dx2-dx1-sdiag1*xp(i-1))*diagin 
!         yp(i) = (dy2-dy1-sdiag1*yp(i-1))*diagin 
!         temp(i) = sdiag2*diagin 
!         dx1 = dx2 
!         dy1 = dy2 
!         diag1 = diag2 
!     8   sdiag1 = sdiag2 
!     9 diag = diag1-sdiag1*temp(nm1) 
!       xp(n) = (slppnx-dx1-sdiag1*xp(nm1))/diag 
!       yp(n) = (slppny-dy1-sdiag1*yp(nm1))/diag 
! !                                                                       
! ! perform back substitution                                             
! !                                                                       
!       do 10 i = 2,n 
!         ibak = np1-i 
!         xp(ibak) = xp(ibak)-temp(ibak)*xp(ibak+1) 
!    10   yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1) 
!       return 
! !                                                                       
! ! too few points                                                        
! !                                                                       
!    11 ierr = 1 
!       return 
! !                                                                       
! ! coincident adjacent points                                            
! !                                                                       
!    12 ierr = 2 
!       return 
!       END SUBROUTINE                                           
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
!       subroutine kurv2 (t,xs,ys,n,x,y,xp,yp,s,sigma) 
! !                                                                       
!       integer n 
!       real t,xs,ys,x(n),y(n),xp(n),yp(n),s(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine performs the mapping of points in the                 
! ! interval (0.,1.) onto a curve in the plane. the subroutine            
! ! kurv1 should be called earlier to determine certain                   
! ! necessary parameters. the resulting curve has a parametric            
! ! representation both of whose components are splines under             
! ! tension and functions of the polygonal arclength                      
! ! parameter.                                                            
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a real value to be mapped to a point on the              
! !   curve. the interval (0.,1.) is mapped onto the entire               
! !   curve, with 0. mapping to (x(1),y(1)) and 1. mapping                
! !   to (x(n),y(n)). values outside this interval result in              
! !   extrapolation.                                                      
! !                                                                       
! !   n contains the number of points which were specified                
! !   to determine the curve.                                             
! !                                                                       
! !   x and y are arrays containing the x- and y-coordinates              
! !   of the specified points.                                            
! !                                                                       
! !   xp and yp are the arrays output from kurv1 containing               
! !   curvature information.                                              
! !                                                                       
! !   s is an array containing the polygonal arclengths of                
! !   the curve.                                                          
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, xp, yp, s, and sigma should be                
! ! input unaltered from the output of kurv1.                             
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   xs and ys contain the x- and y-coordinates of the aimage             
! !   point on the curve.                                                 
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this subroutine references package modules intrvl and                 
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       tn = s(n)*t 
!       im1 = intrvl(tn,s,n) 
!       i = im1+1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/s(n) 
! !                                                                       
! ! set up and perform interpolation                                      
! !                                                                       
!       del1 = tn-s(im1) 
!       del2 = s(i)-tn 
!       dels = s(i)-s(im1) 
!       sumx = (x(i)*del1+x(im1)*del2)/dels 
!       sumy = (y(i)*del1+y(im1)*del2)/dels 
!       if (sigmap .ne. 0.) go to 1 
!       d = del1*del2/(6.*dels) 
!       c1 = (del1+dels)*d 
!       c2 = (del2+dels)*d 
!       xs = sumx-xp(i)*c1-xp(im1)*c2 
!       ys = sumy-yp(i)*c1-yp(im1)*c2 
!       return 
!     1 sigdel = sigmap*dels 
!       call snhcsh(ss,dummy,sigdel,-1) 
!       call snhcsh(s1,dummy,sigmap*del1,-1) 
!       call snhcsh(s2,dummy,sigmap*del2,-1) 
!       d = sigdel*sigmap*(1.+ss) 
!       c1 = del1*(s1-ss)/d 
!       c2 = del2*(s2-ss)/d 
!       xs = sumx+xp(i)*c1+xp(im1)*c2 
!       ys = sumy+yp(i)*c1+yp(im1)*c2 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                         
!       subroutine kurvd (t,xs,ys,xst,yst,xstt,ystt,n,x,y,xp,             &
!      &                  yp,s,sigma)                                     
! !                                                                       
!       integer n 
!       real t,xs,ys,xst,yst,xstt,ystt,x(n),y(n),xp(n),yp(n),             &
!      &     s(n),sigma                                                   
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine performs the mapping of points in the                 
! ! interval (0.,1.) onto a curve in the plane. it also                   
! ! returns the first and second derivatives of the component             
! ! functions. the subroutine kurv1 should be called earlier              
! ! to determine certain necessary parameters. the resulting              
! ! curve has a parametric representation both of whose                   
! ! components are splines under tension and functions of the             
! ! polygonal arclength parameter.                                        
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a real value to be mapped to a point on the              
! !   curve. the interval (0.,1.) is mapped onto the entire               
! !   curve, with 0. mapping to (x(1),y(1)) and 1. mapping                
! !   to (x(n),y(n)). values outside this interval result in              
! !   extrapolation.                                                      
! !                                                                       
! !   n contains the number of points which were specified                
! !   to determine the curve.                                             
! !                                                                       
! !   x and y are arrays containing the x- and y-coordinates              
! !   of the specified points.                                            
! !                                                                       
! !   xp and yp are the arrays output from kurv1 containing               
! !   curvature information.                                              
! !                                                                       
! !   s is an array containing the polygonal arclengths of                
! !   the curve.                                                          
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, xp, yp, s, and sigma should be                
! ! input unaltered from the output of kurv1.                             
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   xs and ys contain the x- and y-coordinates of the aimage             
! !   point on the curve. xst and yst contain the first                   
! !   derivatives of the x- and y-components of the mapping               
! !   with respect to t. xstt and ystt contain the second                 
! !   derivatives of the x- and y-components of the mapping               
! !   with respect to t.                                                  
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this subroutine references package modules intrvl and                 
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       tn = s(n)*t 
!       im1 = intrvl(tn,s,n) 
!       i = im1+1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/s(n) 
! !                                                                       
! ! set up and perform interpolation                                      
! !                                                                       
!       del1 = tn-s(im1) 
!       del2 = s(i)-tn 
!       dels = s(i)-s(im1) 
!       sumx = (x(i)*del1+x(im1)*del2)/dels 
!       sumy = (y(i)*del1+y(im1)*del2)/dels 
!       sumxt = s(n)*(x(i)-x(im1))/dels 
!       sumyt = s(n)*(y(i)-y(im1))/dels 
!       if (sigmap .ne. 0.) go to 1 
!       dels6 = 6.*dels 
!       d = del1*del2/dels6 
!       c1 = -(del1+dels)*d 
!       c2 = -(del2+dels)*d 
!       dels6 = dels6/s(n) 
!       ct1 = (2.*del1*del1-del2*(del1+dels))/dels6 
!       ct2 = -(2.*del2*del2-del1*(del2+dels))/dels6 
!       dels = dels/(s(n)*s(n)) 
!       ctt1 = del1/dels 
!       ctt2 = del2/dels 
!       go to 2 
!     1 sigdel = sigmap*dels 
!       call snhcsh (ss,dummy,sigdel,-1) 
!       call snhcsh (s1,co1,sigmap*del1,0) 
!       call snhcsh (s2,co2,sigmap*del2,0) 
!       d = sigdel*sigmap*(1.+ss) 
!       c1 = del1*(s1-ss)/d 
!       c2 = del2*(s2-ss)/d 
!       ct1 = (co1-ss)*s(n)/d 
!       ct2 = -(co2-ss)*s(n)/d 
!       ctt1 = del1*(1.+s1)*s(n)*s(n)/(dels*(1.+ss)) 
!       ctt2 = del2*(1.+s2)*s(n)*s(n)/(dels*(1.+ss)) 
!     2 xs = sumx+c1*xp(i)+c2*xp(im1) 
!       ys = sumy+c1*yp(i)+c2*yp(im1) 
!       xst = sumxt+ct1*xp(i)+ct2*xp(im1) 
!       yst = sumyt+ct1*yp(i)+ct2*yp(im1) 
!       xstt = ctt1*xp(i)+ctt2*xp(im1) 
!       ystt = ctt1*yp(i)+ctt2*yp(im1) 
!       return 
!       END SUBROUTINE 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                           
!       subroutine kurvp1 (n,x,y,xp,yp,temp,s,sigma,ierr) 
! !                                                                       
!       integer n,ierr 
!       real x(n),y(n),xp(n),yp(n),temp(1),s(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a spline under tension forming a closed curve in              
! ! the plane and passing through a sequence of pairs                     
! ! (x(1),y(1)),...,(x(n),y(n)). for actual computation of                
! ! points on the curve it is necessary to call the subroutine            
! ! kurvp2.                                                               
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of points to be interpolated (n.ge.2).              
! !                                                                       
! !   x is an array containing the n x-coordinates of the                 
! !   points.                                                             
! !                                                                       
! !   y is an array containing the n y-coordinates of the                 
! !   points. (adjacent x-y pairs must be distinct, i. e.                 
! !   either x(i) .ne. x(i+1) or y(i) .ne. y(i+1), for                    
! !   i = 1,...,n-1 and either x(1) .ne. x(n) or y(1) .ne. y(n).)         
! !                                                                       
! !   xp and yp are arrays of length at least n.                          
! !                                                                       
! !   temp is an array of length at least 2*n which is used               
! !   for scratch storage.                                                
! !                                                                       
! !   s is an array of length at least n.                                 
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the curviness desired. if abs(sigma) is nearly zero                 
! !   (e.g. .001) the resulting curve is approximately a cubic            
! !   spline. if abs(sigma) is large (e. g. 50.) the resulting            
! !   curve is nearly a polygonal line. if sigma equals zero a            
! !   cubic spline results. a standard value for sigma is                 
! !   approximately 1. in absolute value.                                 
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   xp and yp contain information about the curvature of the            
! !   curve at the given nodes.                                           
! !                                                                       
! !   s contains the polygonal arclengths of the curve.                   
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if adjacent coordinate pairs coincide.                     
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, and sigma are unaltered,                                   
! !                                                                       
! ! this subroutine references package modules terms and                  
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       nm1 = n-1 
!       np1 = n+1 
!       ierr = 0 
!       if (n .le. 1) go to 7 
! !                                                                       
! ! determine polygonal arclengths                                        
! !                                                                       
!       s(1) = sqrt((x(n)-x(1))**2+(y(n)-y(1))**2) 
!       do 1 i = 2,n 
!         im1 = i-1 
!     1   s(i) = s(im1)+sqrt((x(i)-x(im1))**2+                            &
!      &         (y(i)-y(im1))**2)                                        
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/s(n) 
! !                                                                       
! ! set up right hand sides of tridiagonal (with corner                   
! ! elements) linear system for xp and yp                                 
! !                                                                       
!       dels1 = s(1) 
!       if (dels1 .eq. 0.) go to 8 
!       dx1 = (x(1)-x(n))/dels1 
!       dy1 = (y(1)-y(n))/dels1 
!       call terms(diag1,sdiag1,sigmap,dels1) 
!       dels2 = s(2)-s(1) 
!       if (dels2 .eq. 0.) go to 8 
!       dx2 = (x(2)-x(1))/dels2 
!       dy2 = (y(2)-y(1))/dels2 
!       call terms(diag2,sdiag2,sigmap,dels2) 
!       diag = diag1+diag2 
!       diagin = 1./diag 
!       xp(1) = (dx2-dx1)*diagin 
!       yp(1) = (dy2-dy1)*diagin 
!       temp(np1) = -sdiag1*diagin 
!       temp(1) = sdiag2*diagin 
!       dx1 = dx2 
!       dy1 = dy2 
!       diag1 = diag2 
!       sdiag1 = sdiag2 
!       if (n .eq. 2) go to 3 
!       do 2 i = 2,nm1 
!         npi = n+i 
!         dels2 = s(i+1)-s(i) 
!         if (dels2 .eq. 0.) go to 8 
!         dx2 = (x(i+1)-x(i))/dels2 
!         dy2 = (y(i+1)-y(i))/dels2 
!         call terms(diag2,sdiag2,sigmap,dels2) 
!         diag = diag1+diag2-sdiag1*temp(i-1) 
!         diagin = 1./diag 
!         xp(i) = (dx2-dx1-sdiag1*xp(i-1))*diagin 
!         yp(i) = (dy2-dy1-sdiag1*yp(i-1))*diagin 
!         temp(npi) = -temp(npi-1)*sdiag1*diagin 
!         temp(i) = sdiag2*diagin 
!         dx1 = dx2 
!         dy1 = dy2 
!         diag1 = diag2 
!     2   sdiag1 = sdiag2 
!     3 dels2 = s(1) 
!       dx2 = (x(1)-x(n))/dels2 
!       dy2 = (y(1)-y(n))/dels2 
!       call terms(diag2,sdiag2,sigmap,dels2) 
!       xp(n) = dx2-dx1 
!       yp(n) = dy2-dy1 
!       temp(nm1) = temp(2*n-1)-temp(nm1) 
!       if (n.eq.2) go to 5 
! !                                                                       
! ! perform first step of back substitution                               
! !                                                                       
!       do 4 i = 3,n 
!         ibak = np1-i 
!         npibak = n+ibak 
!         xp(ibak) = xp(ibak)-temp(ibak)*xp(ibak+1) 
!         yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1) 
!     4   temp(ibak) = temp(npibak)-temp(ibak)*temp(ibak+1) 
!     5 xp(n) = (xp(n)-sdiag2*xp(1)-sdiag1*xp(nm1))/                      &
!      &        (diag1+diag2+sdiag2*temp(1)+sdiag1*temp(nm1))             
!       yp(n) = (yp(n)-sdiag2*yp(1)-sdiag1*yp(nm1))/                      &
!      &        (diag1+diag2+sdiag2*temp(1)+sdiag1*temp(nm1))             
! !                                                                       
! ! perform second step of back substitution                              
! !                                                                       
!       xpn = xp(n) 
!       ypn = yp(n) 
!       do 6 i = 1,nm1 
!         xp(i) = xp(i)+temp(i)*xpn 
!     6   yp(i) = yp(i)+temp(i)*ypn 
!       return 
! !                                                                       
! ! too few points                                                        
! !                                                                       
!     7 ierr = 1 
!       return 
! !                                                                       
! ! coincident adjacent points                                            
! !                                                                       
!     8 ierr = 2 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                         
!       subroutine kurvp2 (t,xs,ys,n,x,y,xp,yp,s,sigma) 
! !                                                                       
!       integer n 
!       real t,xs,ys,x(n),y(n),xp(n),yp(n),s(n),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine performs the mapping of points in the                 
! ! interval (0.,1.) onto a closed curve in the plane. the                
! ! subroutine kurvp1 should be called earlier to determine               
! ! certain necessary parameters. the resulting curve has a               
! ! parametric representation both of whose components are                
! ! periodic splines under tension and functions of the poly-             
! ! gonal arclength parameter.                                            
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a value to be mapped onto the curve. the                 
! !   interval (0.,1.) is mapped onto the entire closed curve             
! !   with both 0. and 1. mapping to (x(1),y(1)). the mapping             
! !   is periodic with period one thus any interval of the                
! !   form (tt,tt+1.) maps onto the entire curve.                         
! !                                                                       
! !   n contains the number of points which were specified                
! !   to determine the curve.                                             
! !                                                                       
! !   x and y are arrays containing the x- and y-coordinates              
! !   of the specified points.                                            
! !                                                                       
! !   xp and yp are the arrays output from kurvp1 containing              
! !   curvature information.                                              
! !                                                                       
! !   s is an array containing the polygonal arclengths of                
! !   the curve.                                                          
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, xp, yp, s and sigma should                    
! ! be input unaltered from the output of kurvp1.                         
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   xs and ys contain the x- and y-coordinates of the aimage             
! !   point on the curve.                                                 
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this subroutine references package modules intrvl and                 
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       tn = t-float(ifix(t)) 
!       if (tn .lt. 0.) tn = tn+1. 
!       tn = s(n)*tn+s(1) 
!       im1 = n 
!       if (tn .lt. s(n)) im1 = intrvl(tn,s,n) 
!       i = im1+1 
!       if (i .gt. n) i = 1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/s(n) 
! !                                                                       
! ! set up and perform interpolation                                      
! !                                                                       
!       si = s(i) 
!       if (im1 .eq. n) si = s(n)+s(1) 
!       del1 = tn-s(im1) 
!       del2 = si-tn 
!       dels = si-s(im1) 
!       sumx = (x(i)*del1+x(im1)*del2)/dels 
!       sumy = (y(i)*del1+y(im1)*del2)/dels 
!       if (sigmap .ne. 0.) go to 1 
!       d = del1*del2/(6.*dels) 
!       c1 = (del1+dels)*d 
!       c2 = (del2+dels)*d 
!       xs = sumx-xp(i)*c1-xp(im1)*c2 
!       ys = sumy-yp(i)*c1-yp(im1)*c2 
!       return 
!     1 sigdel = sigmap*dels 
!       call snhcsh(ss,dummy,sigdel,-1) 
!       call snhcsh(s1,dummy,sigmap*del1,-1) 
!       call snhcsh(s2,dummy,sigmap*del2,-1) 
!       d = sigdel*sigmap*(1.+ss) 
!       ci = del1*(s1-ss)/d 
!       cim1 = del2*(s2-ss)/d 
!       xs = sumx+xp(i)*ci+xp(im1)*cim1 
!       ys = sumy+yp(i)*ci+yp(im1)*cim1 
!       return 
!       END SUBROUTINE 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                           
!       subroutine kurvpd (t,xs,ys,xst,yst,xstt,ystt,n,x,y,xp,            &
!      &                  yp,s,sigma)                                     
! !                                                                       
!       integer n 
!       real t,xs,ys,xst,yst,xstt,ystt,x(n),y(n),xp(n),yp(n),             &
!      &     s(n),sigma                                                   
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine performs the mapping of points in the                 
! ! interval (0.,1.) onto a closed curve in the plane. it also            
! ! returns the first and second derivatives of the component             
! ! functions. the subroutine kurvp1 should be called earlier             
! ! to determine certain necessary parameters. the resulting              
! ! curve has a parametric representation both of whose                   
! ! components are periodic splines under tension and                     
! ! functions of the polygonal arclength parameter.                       
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t contains a value to be mapped onto the curve. the                 
! !   interval (0.,1.) is mapped onto the entire closed curve             
! !   with both 0. and 1. mapping to (x(1),y(1)). the mapping             
! !   is periodic with period one thus any interval of the                
! !   form (tt,tt+1.) maps onto the entire curve.                         
! !                                                                       
! !   n contains the number of points which were specified                
! !   to determine the curve.                                             
! !                                                                       
! !   x and y are arrays containing the x- and y-coordinates              
! !   of the specified points.                                            
! !                                                                       
! !   xp and yp are the arrays output from kurvp1 containing              
! !   curvature information.                                              
! !                                                                       
! !   s is an array containing the polygonal arclengths of                
! !   the curve.                                                          
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters n, x, y, xp, yp, s and sigma should                    
! ! be input unaltered from the output of kurvp1.                         
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   xs and ys contain the x- and y-coordinates of the aimage             
! !   point on the curve. xst and yst contain the first                   
! !   derivatives of the x- and y-components of the mapping               
! !   with respect to t. xstt and ystt contain the second                 
! !   derivatives of the x- and y-components of the mapping               
! !   with respect to t.                                                  
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this subroutine references package modules intrvl and                 
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! determine interval                                                    
! !                                                                       
!       tn = t-float(ifix(t)) 
!       if (tn .lt. 0.) tn = tn+1. 
!       tn = s(n)*tn+s(1) 
!       im1 = n 
!       if (tn .lt. s(n)) im1 = intrvl(tn,s,n) 
!       i = im1+1 
!       if (i .gt. n) i = 1 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/s(n) 
! !                                                                       
! ! set up and perform interpolation                                      
! !                                                                       
!       si = s(i) 
!       if (im1 .eq. n) si = s(n)+s(1) 
!       del1 = tn-s(im1) 
!       del2 = si-tn 
!       dels = si-s(im1) 
!       sumx = (x(i)*del1+x(im1)*del2)/dels 
!       sumy = (y(i)*del1+y(im1)*del2)/dels 
!       sumxt = s(n)*(x(i)-x(im1))/dels 
!       sumyt = s(n)*(y(i)-y(im1))/dels 
!       if (sigmap .ne. 0.) go to 1 
!       dels6 = 6.*dels 
!       d = del1*del2/dels6 
!       c1 = -(del1+dels)*d 
!       c2 = -(del2+dels)*d 
!       dels6 = dels6/s(n) 
!       ct1 = (2.*del1*del1-del2*(del1+dels))/dels6 
!       ct2 = -(2.*del2*del2-del1*(del2+dels))/dels6 
!       dels = dels/(s(n)*s(n)) 
!       ctt1 = del1/dels 
!       ctt2 = del2/dels 
!       go to 2 
!     1 sigdel = sigmap*dels 
!       call snhcsh (ss,dummy,sigdel,-1) 
!       call snhcsh (s1,co1,sigmap*del1,0) 
!       call snhcsh (s2,co2,sigmap*del2,0) 
!       d = sigdel*sigmap*(1.+ss) 
!       c1 = del1*(s1-ss)/d 
!       c2 = del2*(s2-ss)/d 
!       ct1 = (co1-ss)*s(n)/d 
!       ct2 = -(co2-ss)*s(n)/d 
!       ctt1 = del1*(1.+s1)*s(n)*s(n)/(dels*(1.+ss)) 
!       ctt2 = del2*(1.+s2)*s(n)*s(n)/(dels*(1.+ss)) 
!     2 xs = sumx+c1*xp(i)+c2*xp(im1) 
!       ys = sumy+c1*yp(i)+c2*yp(im1) 
!       xst = sumxt+ct1*xp(i)+ct2*xp(im1) 
!       yst = sumyt+ct1*yp(i)+ct2*yp(im1) 
!       xstt = ctt1*xp(i)+ctt2*xp(im1) 
!       ystt = ctt1*yp(i)+ctt2*yp(im1) 
!       return 
!       END SUBROUTINE
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!                                            
!       subroutine surf1 (m,n,x,y,z,iz,zx1,zxm,zy1,zyn,zxy11,             &
!      &                  zxym1,zxy1n,zxymn,islpsw,zp,temp,               &
!      &                  sigma,ierr)                                     
! !                                                                       
!       integer m,n,iz,islpsw,ierr 
!       real x(m),y(n),z(iz,n),zx1(n),zxm(n),zy1(m),zyn(m),               &
!      &     zxy11,zxym1,zxy1n,zxymn,zp(m,n,3),temp(1),sigma              
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute an interpolatory surface passing through a rect-              
! ! angular grid of functional values. the surface determined             
! ! can be represented as the tensor product of splines under             
! ! tension. the x- and y-partial derivatives around the                  
! ! boundary and the x-y-partial derivatives at the four                  
! ! corners may be specified or omitted. for actual mapping               
! ! of points onto the surface it is necessary to call the                
! ! function surf2.                                                       
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   m is the number of grid lines in the x-direction, i. e.             
! !   lines parallel to the y-axis (m .ge. 2).                            
! !                                                                       
! !   n is the number of grid lines in the y-direction, i. e.             
! !   lines parallel to the x-axis (n .ge. 2).                            
! !                                                                       
! !   x is an array of the m x-coordinates of the grid lines              
! !   in the x-direction. these should be strictly increasing.            
! !                                                                       
! !   y is an array of the n y-coordinates of the grid lines              
! !   in the y-direction. these should be strictly increasing.            
! !                                                                       
! !   z is an array of the m * n functional values at the grid            
! !   points, i. e. z(i,j) contains the functional value at               
! !   (x(i),y(j)) for i = 1,...,m and j = 1,...,n.                        
! !                                                                       
! !   iz is the row dimension of the matrix z used in the                 
! !   calling program (iz .ge. m).                                        
! !                                                                       
! !   zx1 and zxm are arrays of the m x-partial derivatives               
! !   of the function along the x(1) and x(m) grid lines,                 
! !   respectively. thus zx1(j) and zxm(j) contain the x-part-            
! !   ial derivatives at the points (x(1),y(j)) and                       
! !   (x(m),y(j)), respectively, for j = 1,...,n. either of               
! !   these parameters will be ignored (and approximations                
! !   supplied internally) if islpsw so indicates.                        
! !                                                                       
! !   zy1 and zyn are arrays of the n y-partial derivatives               
! !   of the function along the y(1) and y(n) grid lines,                 
! !   respectively. thus zy1(i) and zyn(i) contain the y-part-            
! !   ial derivatives at the points (x(i),y(1)) and                       
! !   (x(i),y(n)), respectively, for i = 1,...,m. either of               
! !   these parameters will be ignored (and estimations                   
! !   supplied internally) if islpsw so indicates.                        
! !                                                                       
! !   zxy11, zxym1, zxy1n, and zxymn are the x-y-partial                  
! !   derivatives of the function at the four corners,                    
! !   (x(1),y(1)), (x(m),y(1)), (x(1),y(n)), and (x(m),y(n)),             
! !   respectively. any of the parameters will be ignored (and            
! !   estimations supplied internally) if islpsw so indicates.            
! !                                                                       
! !   islpsw contains a switch indicating which boundary                  
! !   derivative information is user-supplied and which                   
! !   should be estimated by this subroutine. to determine                
! !   islpsw, let                                                         
! !        i1 = 0 if zx1 is user-supplied (and = 1 otherwise),            
! !        i2 = 0 if zxm is user-supplied (and = 1 otherwise),            
! !        i3 = 0 if zy1 is user-supplied (and = 1 otherwise),            
! !        i4 = 0 if zyn is user-supplied (and = 1 otherwise),            
! !        i5 = 0 if zxy11 is user-supplied                               
! !                                       (and = 1 otherwise),            
! !        i6 = 0 if zxym1 is user-supplied                               
! !                                       (and = 1 otherwise),            
! !        i7 = 0 if zxy1n is user-supplied                               
! !                                       (and = 1 otherwise),            
! !        i8 = 0 if zxymn is user-supplied                               
! !                                       (and = 1 otherwise),            
! !   then islpsw = i1 + 2*i2 + 4*i3 + 8*i4 + 16*i5 + 32*i6               
! !                   + 64*i7 + 128*i8                                    
! !   thus islpsw = 0 indicates all derivative information is             
! !   user-supplied and islpsw = 255 indicates no derivative              
! !   information is user-supplied. any value between these               
! !   limits is valid.                                                    
! !                                                                       
! !   zp is an array of at least 3*m*n locations.                         
! !                                                                       
! !   temp is an array of at least n+n+m locations which is               
! !   used for scratch storage.                                           
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the curviness desired. if abs(sigma) is nearly zero                 
! !   (e. g. .001) the resulting surface is approximately the             
! !   tensor product of cubic splines. if abs(sigma) is large             
! !   (e. g. 50.) the resulting surface is approximately                  
! !   bi-linear. if sigma equals zero tensor products of                  
! !   cubic splines result. a standard value for sigma is                 
! !   approximately 1. in absolute value.                                 
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   zp contains the values of the xx-, yy-, and xxyy-partial            
! !   derivatives of the surface at the given nodes.                      
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2 or m is less than 2,                   
! !        = 2 if the x-values or y-values are not strictly               
! !            increasing.                                                
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1,                
! !   zxy1n, zxymn, islpsw, and sigma are unaltered.                      
! !                                                                       
! ! this subroutine references package modules ceez, terms,               
! ! and snhcsh.                                                           
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       mm1 = m-1 
!       mp1 = m+1 
!       nm1 = n-1 
!       np1 = n+1 
!       npm = n+m 
!       ierr = 0 
!       if (n .le. 1 .or. m .le. 1) go to 46 
!       if (y(n) .le. y(1)) go to 47 
! !                                                                       
! ! denormalize tension factor in y-direction                             
! !                                                                       
!       sigmay = abs(sigma)*float(n-1)/(y(n)-y(1)) 
! !                                                                       
! ! obtain y-partial derivatives along y = y(1)                           
! !                                                                       
!       if ((islpsw/8)*2 .ne. (islpsw/4)) go to 2 
!       do 1 i = 1,m 
!     1   zp(i,1,1) = zy1(i) 
!       go to 5 
!     2 dely1 = y(2)-y(1) 
!       dely2 = dely1+dely1 
!       if (n .gt. 2) dely2 = y(3)-y(1) 
!       if (dely1 .le. 0. .or. dely2 .le. dely1) go to 47 
!       call ceez (dely1,dely2,sigmay,c1,c2,c3,n) 
!       do 3 i = 1,m 
!     3   zp(i,1,1) = c1*z(i,1)+c2*z(i,2) 
!       if (n .eq. 2) go to 5 
!       do 4 i = 1,m 
!     4   zp(i,1,1) = zp(i,1,1)+c3*z(i,3) 
! !                                                                       
! ! obtain y-partial derivatives along y = y(n)                           
! !                                                                       
!     5 if ((islpsw/16)*2 .ne. (islpsw/8)) go to 7 
!       do 6 i = 1,m 
!         npi = n+i 
!     6   temp(npi) = zyn(i) 
!       go to 10 
!     7 delyn = y(n)-y(nm1) 
!       delynm = delyn+delyn 
!       if (n .gt. 2) delynm = y(n)-y(n-2) 
!       if (delyn .le. 0. .or. delynm .le. delyn) go to 47 
!       call ceez (-delyn,-delynm,sigmay,c1,c2,c3,n) 
!       do 8 i = 1,m 
!         npi = n+i 
!     8   temp(npi) = c1*z(i,n)+c2*z(i,nm1) 
!       if (n .eq. 2) go to 10 
!       do 9 i = 1,m 
!         npi = n+i 
!     9   temp(npi) = temp(npi)+c3*z(i,n-2) 
!    10 if (x(m) .le. x(1)) go to 47 
! !                                                                       
! ! denormalize tension factor in x-direction                             
! !                                                                       
!       sigmax = abs(sigma)*float(m-1)/(x(m)-x(1)) 
! !                                                                       
! ! obtain x-partial derivatives along x = x(1)                           
! !                                                                       
!       if ((islpsw/2)*2 .ne. islpsw) go to 12 
!       do 11 j = 1,n 
!    11   zp(1,j,2) = zx1(j) 
!       if ((islpsw/32)*2 .eq. (islpsw/16) .and.                          &
!      &    (islpsw/128)*2  .eq. (islpsw/64)) go to 15                    
!    12 delx1 = x(2)-x(1) 
!       delx2 = delx1+delx1 
!       if (m .gt. 2) delx2 = x(3)-x(1) 
!       if (delx1 .le. 0. .or. delx2 .le. delx1) go to 47 
!       call ceez (delx1,delx2,sigmax,c1,c2,c3,m) 
!       if ((islpsw/2)*2 .eq. islpsw) go to 15 
!       do 13 j = 1,n 
!    13   zp(1,j,2) = c1*z(1,j)+c2*z(2,j) 
!       if (m .eq. 2) go to 15 
!       do 14 j = 1,n 
!    14   zp(1,j,2) = zp(1,j,2)+c3*z(3,j) 
! !                                                                       
! ! obtain x-y-partial derivative at (x(1),y(1))                          
! !                                                                       
!    15 if ((islpsw/32)*2 .ne. (islpsw/16)) go to 16 
!       zp(1,1,3) = zxy11 
!       go to 17 
!    16 zp(1,1,3) = c1*zp(1,1,1)+c2*zp(2,1,1) 
!       if (m .gt. 2) zp(1,1,3) = zp(1,1,3)+c3*zp(3,1,1) 
! !                                                                       
! ! obtain x-y-partial derivative at (x(1),y(n))                          
! !                                                                       
!    17 if ((islpsw/128)*2 .ne. (islpsw/64)) go to 18 
!       zxy1ns = zxy1n 
!       go to 19 
!    18 zxy1ns = c1*temp(n+1)+c2*temp(n+2) 
!       if (m .gt. 2) zxy1ns = zxy1ns+c3*temp(n+3) 
! !                                                                       
! ! obtain x-partial derivative along x = x(m)                            
! !                                                                       
!    19 if ((islpsw/4)*2 .ne. (islpsw/2)) go to 21 
!       do 20 j = 1,n 
!         npmpj = npm+j 
!    20   temp(npmpj) = zxm(j) 
!       if ((islpsw/64)*2 .eq. (islpsw/32) .and.                          &
!      &    (islpsw/256)*2 .eq. (islpsw/128)) go to 24                    
!    21 delxm = x(m)-x(mm1) 
!       delxmm = delxm+delxm 
!       if (m .gt. 2) delxmm = x(m)-x(m-2) 
!       if (delxm .le. 0. .or. delxmm .le. delxm) go to 47 
!       call ceez (-delxm,-delxmm,sigmax,c1,c2,c3,m) 
!       if ((islpsw/4)*2 .eq. (islpsw/2)) go to 24 
!       do 22 j = 1,n 
!         npmpj = npm+j 
!    22   temp(npmpj) = c1*z(m,j)+c2*z(mm1,j) 
!       if (m .eq. 2) go to 24 
!       do 23 j = 1,n 
!         npmpj = npm+j 
!    23   temp(npmpj) = temp(npmpj)+c3*z(m-2,j) 
! !                                                                       
! ! obtain x-y-partial derivative at (x(m),y(1))                          
! !                                                                       
!    24 if ((islpsw/64)*2 .ne. (islpsw/32)) go to 25 
!       zp(m,1,3) = zxym1 
!       go to 26 
!    25 zp(m,1,3) = c1*zp(m,1,1)+c2*zp(mm1,1,1) 
!       if (m .gt. 2) zp(m,1,3) = zp(m,1,3)+c3*zp(m-2,1,1) 
! !                                                                       
! ! obtain x-y-partial derivative at (x(m),y(n))                          
! !                                                                       
!    26 if ((islpsw/256)*2 .ne. (islpsw/128)) go to 27 
!       zxymns = zxymn 
!       go to 28 
!    27 zxymns = c1*temp(npm)+c2*temp(npm-1) 
!       if (m .gt. 2) zxymns = zxymns+c3*temp(npm-2) 
! !                                                                       
! ! set up right hand sides and tridiagonal system for y-grid             
! ! perform forward elimination                                           
! !                                                                       
!    28 del1 = y(2)-y(1) 
!       if (del1 .le. 0.) go to 47 
!       deli = 1./del1 
!       do 29 i = 1,m 
!    29   zp(i,2,1) = deli*(z(i,2)-z(i,1)) 
!       zp(1,2,3) = deli*(zp(1,2,2)-zp(1,1,2)) 
!       zp(m,2,3) = deli*(temp(npm+2)-temp(npm+1)) 
!       call terms (diag1,sdiag1,sigmay,del1) 
!       diagi = 1./diag1 
!       do 30 i = 1,m 
!    30   zp(i,1,1) = diagi*(zp(i,2,1)-zp(i,1,1)) 
!       zp(1,1,3) = diagi*(zp(1,2,3)-zp(1,1,3)) 
!       zp(m,1,3) = diagi*(zp(m,2,3)-zp(m,1,3)) 
!       temp(1) = diagi*sdiag1 
!       if (n .eq. 2) go to 34 
!       do 33 j = 2,nm1 
!         jm1 = j-1 
!         jp1 = j+1 
!         npmpj = npm+j 
!         del2 = y(jp1)-y(j) 
!         if (del2 .le. 0.) go to 47 
!         deli = 1./del2 
!         do 31 i = 1,m 
!    31     zp(i,jp1,1) = deli*(z(i,jp1)-z(i,j)) 
!         zp(1,jp1,3) = deli*(zp(1,jp1,2)-zp(1,j,2)) 
!         zp(m,jp1,3) = deli*(temp(npmpj+1)-temp(npmpj)) 
!         call terms (diag2,sdiag2,sigmay,del2) 
!         diagin = 1./(diag1+diag2-sdiag1*temp(jm1)) 
!         do 32 i = 1,m 
!    32     zp(i,j,1) = diagin*(zp(i,jp1,1)-zp(i,j,1)-                    &
!      &                        sdiag1*zp(i,jm1,1))                       
!         zp(1,j,3) = diagin*(zp(1,jp1,3)-zp(1,j,3)-                      &
!      &                      sdiag1*zp(1,jm1,3))                         
!         zp(m,j,3) = diagin*(zp(m,jp1,3)-zp(m,j,3)-                      &
!      &                      sdiag1*zp(m,jm1,3))                         
!         temp(j) = diagin*sdiag2 
!         diag1 = diag2 
!    33   sdiag1 = sdiag2 
!    34 diagin = 1./(diag1-sdiag1*temp(nm1)) 
!       do 35 i = 1,m 
!         npi = n+i 
!    35   zp(i,n,1) = diagin*(temp(npi)-zp(i,n,1)-                        &
!      &                      sdiag1*zp(i,nm1,1))                         
!       zp(1,n,3) = diagin*(zxy1ns-zp(1,n,3)-                             &
!      &                    sdiag1*zp(1,nm1,3))                           
!       temp(n) = diagin*(zxymns-zp(m,n,3)-                               &
!      &                  sdiag1*zp(m,nm1,3))                             
! !                                                                       
! ! perform back substitution                                             
! !                                                                       
!       do 37 j = 2,n 
!         jbak = np1-j 
!         jbakp1 = jbak+1 
!         t = temp(jbak) 
!         do 36 i = 1,m 
!    36     zp(i,jbak,1) = zp(i,jbak,1)-t*zp(i,jbakp1,1) 
!         zp(1,jbak,3) = zp(1,jbak,3)-t*zp(1,jbakp1,3) 
!    37   temp(jbak) = zp(m,jbak,3)-t*temp(jbakp1) 
! !                                                                       
! ! set up right hand sides and tridiagonal system for x-grid             
! ! perform forward elimination                                           
! !                                                                       
!       del1 = x(2)-x(1) 
!       if (del1 .le. 0.) go to 47 
!       deli = 1./del1 
!       do 38 j = 1,n 
!         zp(2,j,2) = deli*(z(2,j)-z(1,j)) 
!    38   zp(2,j,3) = deli*(zp(2,j,1)-zp(1,j,1)) 
!       call terms (diag1,sdiag1,sigmax,del1) 
!       diagi = 1./diag1 
!       do 39 j = 1,n 
!         zp(1,j,2) = diagi*(zp(2,j,2)-zp(1,j,2)) 
!    39   zp(1,j,3) = diagi*(zp(2,j,3)-zp(1,j,3)) 
!       temp(n+1) = diagi*sdiag1 
!       if (m  .eq. 2) go to 43 
!       do 42 i = 2,mm1 
!         im1 = i-1 
!         ip1 = i+1 
!         npi = n+i 
!         del2 = x(ip1)-x(i) 
!         if (del2 .le. 0.) go to 47 
!         deli = 1./del2 
!         do 40 j = 1,n 
!           zp(ip1,j,2) = deli*(z(ip1,j)-z(i,j)) 
!    40     zp(ip1,j,3) = deli*(zp(ip1,j,1)-zp(i,j,1)) 
!         call terms (diag2,sdiag2,sigmax,del2) 
!         diagin = 1./(diag1+diag2-sdiag1*temp(npi-1)) 
!         do 41 j = 1,n 
!           zp(i,j,2) = diagin*(zp(ip1,j,2)-zp(i,j,2)-                    &
!      &                        sdiag1*zp(im1,j,2))                       
!    41     zp(i,j,3) = diagin*(zp(ip1,j,3)-zp(i,j,3)-                    &
!      &                        sdiag1*zp(im1,j,3))                       
!         temp(npi) = diagin*sdiag2 
!         diag1 = diag2 
!    42   sdiag1 = sdiag2 
!    43 diagin = 1./(diag1-sdiag1*temp(npm-1)) 
!       do 44 j = 1,n 
!         npmpj = npm+j 
!         zp(m,j,2) = diagin*(temp(npmpj)-zp(m,j,2)-                      &
!      &                      sdiag1*zp(mm1,j,2))                         
!    44   zp(m,j,3) = diagin*(temp(j)-zp(m,j,3)-                          &
!      &                      sdiag1*zp(mm1,j,3))                         
! !                                                                       
! ! perform back substitution                                             
! !                                                                       
!       do 45 i = 2,m 
!         ibak = mp1-i 
!         ibakp1 = ibak+1 
!         npibak = n+ibak 
!         t = temp(npibak) 
!         do 45 j = 1,n 
!           zp(ibak,j,2) = zp(ibak,j,2)-t*zp(ibakp1,j,2) 
!    45     zp(ibak,j,3) = zp(ibak,j,3)-t*zp(ibakp1,j,3) 
!       return 
! !                                                                       
! ! too few points                                                        
! !                                                                       
!    46 ierr = 1 
!       return 
! !                                                                       
! ! points not strictly increasing                                        
! !                                                                       
!    47 ierr = 2 
!       return 
!       END SUBROUTINE  
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!                                          
!       function surf2 (xx,yy,m,n,x,y,z,iz,zp,sigma) 
! !                                                                       
!       integer m,n,iz 
!       real xx,yy,x(m),y(n),z(iz,n),zp(m,n,3),sigma 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function interpolates a surface at a given coordinate            
! ! pair using a bi-spline under tension. the subroutine surf1            
! ! should be called earlier to determine certain necessary               
! ! parameters.                                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   xx and yy contain the x- and y-coordinates of the point             
! !   to be mapped onto the interpolating surface.                        
! !                                                                       
! !   m and n contain the number of grid lines in the x- and              
! !   y-directions, respectively, of the rectangular grid                 
! !   which specified the surface.                                        
! !                                                                       
! !   x and y are arrays containing the x- and y-grid values,             
! !   respectively, each in increasing order.                             
! !                                                                       
! !   z is a matrix containing the m * n functional values                
! !   corresponding to the grid values (i. e. z(i,j) is the               
! !   surface value at the point (x(i),y(j)) for i = 1,...,m              
! !   and j = 1,...,n).                                                   
! !                                                                       
! !   iz contains the row dimension of the array z as declared            
! !   in the calling program.                                             
! !                                                                       
! !   zp is an array of 3*m*n locations stored with the                   
! !   various surface derivative information determined by                
! !   surf1.                                                              
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma contains the tension factor (its sign is ignored).            
! !                                                                       
! ! the parameters m, n, x, y, z, iz, zp, and sigma should be             
! ! input unaltered from the output of surf1.                             
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   surf2 contains the interpolated surface value.                      
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this function references package modules intrvl and                   
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
! ! inline one dimensional cubic spline interpolation                     
! !                                                                       
!       hermz (f1,f2,fp1,fp2) = (f2*del1+f1*del2)/dels-del1*              &
!      &                        del2*(fp2*(del1+dels)+                    &
!      &                              fp1*(del2+dels))/                   &
!      &                        (6.*dels)                                 
! !                                                                       
! ! inline one dimensional spline under tension interpolation             
! !                                                                       
!       hermnz (f1,f2,fp1,fp2,sigmap) = (f2*del1+f1*del2)/dels            &
!      &          +(fp2*del1*(sinhm1-sinhms)                              &
!      &           +fp1*del2*(sinhm2-sinhms)                              &
!      &          )/(sigmap*sigmap*dels*(1.+sinhms))                      
! !                                                                       
! ! denormalize tension factor in x and y direction                       
! !                                                                       
!       sigmax = abs(sigma)*float(m-1)/(x(m)-x(1)) 
!       sigmay = abs(sigma)*float(n-1)/(y(n)-y(1)) 
! !                                                                       
! ! determine y interval                                                  
! !                                                                       
!       jm1 = intrvl (yy,y,n) 
!       j = jm1+1 
! !                                                                       
! ! determine x interval                                                  
! !                                                                       
!       im1 = intrvl (xx,x,m) 
!       i = im1+1 
!       del1 = yy-y(jm1) 
!       del2 = y(j)-yy 
!       dels = y(j)-y(jm1) 
!       if (sigmay .ne. 0.) go to 1 
! !                                                                       
! ! perform four interpolations in y-direction                            
! !                                                                       
!       zim1 = hermz(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),                   &
!      &                                  zp(i-1,j,1))                    
!       zi = hermz(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1)) 
!       zxxim1 = hermz(zp(i-1,j-1,2),zp(i-1,j,2),                         &
!      &                zp(i-1,j-1,3),zp(i-1,j,3))                        
!       zxxi = hermz(zp(i,j-1,2),zp(i,j,2),                               &
!      &              zp(i,j-1,3),zp(i,j,3))                              
!       go to 2 
!     1 call snhcsh (sinhm1,dummy,sigmay*del1,-1) 
!       call snhcsh (sinhm2,dummy,sigmay*del2,-1) 
!       call snhcsh (sinhms,dummy,sigmay*dels,-1) 
!       zim1 = hermnz(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),                  &
!      &               zp(i-1,j,1),sigmay)                                
!       zi = hermnz(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1),                &
!      &             sigmay)                                              
!       zxxim1 = hermnz(zp(i-1,j-1,2),zp(i-1,j,2),                        &
!      &                 zp(i-1,j-1,3),zp(i-1,j,3),sigmay)                
!       zxxi = hermnz(zp(i,j-1,2),zp(i,j,2),                              &
!      &               zp(i,j-1,3),zp(i,j,3),sigmay)                      
! !                                                                       
! ! perform final interpolation in x-direction                            
! !                                                                       
!     2 del1 = xx-x(im1) 
!       del2 = x(i)-xx 
!       dels = x(i)-x(im1) 
!       if (sigmax .ne. 0.) go to 3 
!       surf2 = hermz(zim1,zi,zxxim1,zxxi) 
!       return 
!     3 call snhcsh (sinhm1,dummy,sigmax*del1,-1) 
!       call snhcsh (sinhm2,dummy,sigmax*del2,-1) 
!       call snhcsh (sinhms,dummy,sigmax*dels,-1) 
!       surf2 = hermnz(zim1,zi,zxxim1,zxxi,sigmax) 
!       return 
!       END function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!                                         
!       subroutine ceez (del1,del2,sigma,c1,c2,c3,n) 
! !                                                                       
!       real del1,del2,sigma,c1,c2,c3 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the coefficients c1, c2, and c3            
! ! used to determine endpoint slopes. specifically, if                   
! ! function values y1, y2, and y3 are given at points x1, x2,            
! ! and x3, respectively, the quantity c1*y1 + c2*y2 + c3*y3              
! ! is the value of the derivative at x1 of a spline under                
! ! tension (with tension factor sigma) passing through the               
! ! three points and having third derivative equal to zero at             
! ! x1. optionally, only two values, c1 and c2 are determined.            
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   del1 is x2-x1 (.gt. 0.).                                            
! !                                                                       
! !   del2 is x3-x1 (.gt. 0.). if n .eq. 2, this parameter is             
! !   ignored.                                                            
! !                                                                       
! !   sigma is the tension factor.                                        
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n is a switch indicating the number of coefficients to              
! !   be returned. if n .eq. 2 only two coefficients are                  
! !   returned. otherwise all three are returned.                         
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   c1, c2, and c3 contain the coefficients.                            
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! ! this subroutine references package module snhcsh.                     
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       if (n .eq. 2) go to 2 
!       if (sigma .ne. 0.) go to 1 
!       del = del2-del1 
! !                                                                       
! ! tension .eq. 0.                                                       
! !                                                                       
!       c1 = -(del1+del2)/(del1*del2) 
!       c2 = del2/(del1*del) 
!       c3 = -del1/(del2*del) 
!       return 
! !                                                                       
! ! tension .ne. 0.                                                       
! !                                                                       
!     1 call snhcsh (dummy,coshm1,sigma*del1,1) 
!       call snhcsh (dummy,coshm2,sigma*del2,1) 
!       delp = sigma*(del2+del1)/2. 
!       delm = sigma*(del2-del1)/2. 
!       call snhcsh (sinhmp,dummy,delp,-1) 
!       call snhcsh (sinhmm,dummy,delm,-1) 
!       denom = coshm1*(del2-del1)-2.*del1*delp*delm*                     &
!      &        (1.+sinhmp)*(1.+sinhmm)                                   
!       c1 = 2.*delp*delm*(1.+sinhmp)*(1.+sinhmm)/denom 
!       c2 = -coshm2/denom 
!       c3 = coshm1/denom 
!       return 
! !                                                                       
! ! two coefficients                                                      
! !                                                                       
!     2 c1 = -1./del1 
!       c2 = -c1 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                         
!       subroutine curvpp (n,x,y,p,d,isw,s,eps,ys,ysp,sigma,              &
!      &                   td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,             &
!      &                   rnm1,rn,v,ierr)                                
! !                                                                       
!       integer n,isw,ierr 
!       real x(n),y(n),p,d(n),s,eps,ys(n),ysp(n),sigma,td(n),             &
!      &     tsd1(n),hd(n),hsd1(n),hsd2(n),rd(n),rsd1(n),                 &
!      &     rsd2(n),rnm1(n),rn(n),v(n)                                   
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a periodic smoothing spline under tension. for a              
! ! given increasing sequence of abscissae (x(i)), i = 1,...,n            
! ! and associated ordinates (y(i)), i = 1,...,n, letting p be            
! ! the period, x(n+1) = x(1)+p, and y(n+1) = y(1), the                   
! ! function determined minimizes the summation from i = 1 to             
! ! n of the square of the second derivative of f plus sigma              
! ! squared times the difference of the first derivative of f             
! ! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all               
! ! functions f with period p and two continuous derivatives              
! ! such that the summation of the square of                              
! ! (f(x(i))-y(i))/d(i) is less than or equal to a given                  
! ! constant s, where (d(i)), i = 1,...,n are a given set of              
! ! observation weights. the function determined is a periodic            
! ! spline under tension with third derivative discontinuities            
! ! at (x(i)) i = 1,...,n (and all periodic translations of               
! ! these values). for actual computation of points on the                
! ! curve it is necessary to call the function curvp2.                    
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of values to be smoothed (n.ge.2).                  
! !                                                                       
! !   x is an array of the n increasing abscissae of the                  
! !   values to be smoothed.                                              
! !                                                                       
! !   y is an array of the n ordinates of the values to be                
! !   smoothed, (i. e. y(k) is the functional value                       
! !   corresponding to x(k) ).                                            
! !                                                                       
! !   p is the period (p .gt. x(n)-x(1)).                                 
! !                                                                       
! !   d is a parameter containing the observation weights.                
! !   this may either be an array of length n or a scalar                 
! !   (interpreted as a constant). the value of d                         
! !   corresponding to the observation (x(k),y(k)) should                 
! !   be an approximation to the standard deviation of error.             
! !                                                                       
! !   isw contains a switch indicating whether the parameter              
! !   d is to be considered a vector or a scalar,                         
! !          = 0 if d is an array of length n,                            
! !          = 1 if d is a scalar.                                        
! !                                                                       
! !   s contains the value controlling the smoothing. this                
! !   must be non-negative. for s equal to zero, the                      
! !   subroutine does interpolation, larger values lead to                
! !   smoother funtions. if parameter d contains standard                 
! !   deviation estimates, a reasonable value for s is                    
! !   float(n).                                                           
! !                                                                       
! !   eps contains a tolerance on the relative precision to               
! !   which s is to be interpreted. this must be greater than             
! !   or equal to zero and less than equal or equal to one. a             
! !   reasonable value for eps is sqrt(2./float(n)).                      
! !                                                                       
! !   ys is an array of length at least n.                                
! !                                                                       
! !   ysp is an array of length at least n.                               
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the degree to which the first derivative part of the                
! !   smoothing functional is emphasized. if sigma is nearly              
! !   zero (e. g. .001) the resulting curve is approximately a            
! !   cubic spline. if sigma is large (e. g. 50.) the                     
! !   resulting curve is nearly a polygonal line. if sigma                
! !   equals zero a cubic spline results. a standard value for            
! !   sigma is approximately 1.                                           
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, rnm1, rn, and             
! !   v are arrays of length at least n which are used for                
! !   scratch storage.                                                    
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   ys contains the smoothed ordinate values.                           
! !                                                                       
! !   ysp contains the values of the second derivative of the             
! !   smoothed curve at the given nodes.                                  
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if s is negative,                                          
! !        = 3 if eps is negative or greater than one,                    
! !        = 4 if x-values are not strictly increasing,                   
! !        = 5 if a d-value is non-positive,                              
! !        = 6 if p is less than or equal to x(n)-x(1).                   
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, d, isw, s, eps, and sigma are unaltered.                   
! !                                                                       
! ! this subroutine references package modules terms and                  
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       if (n .lt. 2) go to 25 
!       if (s .lt. 0.) go to 26 
!       if (eps .lt. 0. .or. eps .gt. 1.) go to 27 
!       if (p .le. x(n)-x(1)) go to 30 
!       ierr = 0 
!       q = 0. 
!       rsd1(1) = 0. 
!       rsd2(1) = 0. 
!       rsd2(2) = 0. 
!       rsd1(n-1) = 0. 
!       rsd2(n-1) = 0. 
!       rsd2(n) = 0. 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n)/p 
! !                                                                       
! ! form t matrix and second differences of y into ys                     
! !                                                                       
!       nm1 = n-1 
!       nm2 = n-2 
!       nm3 = n-3 
!       delxi1 = x(1)+p-x(n) 
!       delyi1 = (y(1)-y(n))/delxi1 
!       call terms (dim1,tsd1(1),sigmap,delxi1) 
!       hsd1(1) = 1./delxi1 
!       do 1 i = 1,n 
!         ip1 = i+1 
!         if (i .eq. n) ip1 = 1 
!         delxi = x(ip1)-x(i) 
!         if (i .eq. n) delxi = x(1)+p-x(n) 
!         if (delxi .le. 0.) go to 28 
!         delyi = (y(ip1)-y(i))/delxi 
!         ys(i) = delyi-delyi1 
!         call terms (di,tsd1(ip1),sigmap,delxi) 
!         td(i) = di+dim1 
!         hd(i) = -(1./delxi+1./delxi1) 
!         hsd1(ip1) = 1./delxi 
!         delxi1 = delxi 
!         delyi1 = delyi 
!     1   dim1 = di 
!       hsd11 = hsd1(1) 
!       if (n .ge. 3) go to 2 
!       tsd1(2) = tsd1(1)+tsd1(2) 
!       tsd1(1) = 0. 
!       hsd1(2) = hsd1(1)+hsd1(2) 
!       hsd1(1) = 0. 
! !                                                                       
! ! calculate lower and upper tolerances                                  
! !                                                                       
!     2 sl = s*(1.-eps) 
!       su = s*(1.+eps) 
!       if (d(1) .le. 0.) go to 29 
!       if (isw .eq. 1) go to 5 
! !                                                                       
! ! form h matrix - d array                                               
! !                                                                       
!       betapp = hsd1(n)*d(n)*d(n) 
!       betap = hsd1(1)*d(1)*d(1) 
!       alphap = hd(n)*d(n)*d(n) 
!       im1 = n 
!       sumd = 0. 
!       sumy = 0. 
!       do 3 i = 1,n 
!         disq = d(i)*d(i) 
!         sumd = sumd+1./disq 
!         sumy = sumy+y(i)/disq 
!         ip1 = i+1 
!         if (i .eq. n) ip1 = 1 
!         alpha = hd(i)*disq 
!         if (d(ip1) .le. 0.) go to 29 
!         hsd1ip = hsd1(ip1) 
!         if (i .eq. n) hsd1ip = hsd11 
!         beta = hsd1ip*d(ip1)*d(ip1) 
!         hd(i) = (hsd1(i)*d(im1))**2+alpha*hd(i)                         &
!      &                             +beta*hsd1ip                         
!         hsd2(i) = hsd1(i)*betapp 
!         hsd1(i) = hsd1(i)*(alpha+alphap) 
!         im1 = i 
!         alphap = alpha 
!         betapp = betap 
!     3   betap = beta 
!       if (n .eq. 3) hsd1(3) = hsd1(3)+hsd2(2) 
! !                                                                       
! ! test for straight line fit                                            
! !                                                                       
!       con = sumy/sumd 
!       sum = 0. 
!       do 4 i = 1,n 
!     4   sum = sum+((y(i)-con)/d(i))**2 
!       if (sum .le. su) go to 23 
!       go to 8 
! !                                                                       
! ! form h matrix - d constant                                            
! !                                                                       
!     5 sl = d(1)*d(1)*sl 
!       su = d(1)*d(1)*su 
!       hsd1p = hsd1(n) 
!       hdim1 = hd(n) 
!       sumy = 0. 
!       do 6 i = 1,n 
!         sumy = sumy+y(i) 
!         hsd1ip = hsd11 
!         if (i .lt. n) hsd1ip = hsd1(i+1) 
!         hdi = hd(i) 
!         hd(i) = hsd1(i)*hsd1(i)+hdi*hdi+hsd1ip*hsd1ip 
!         hsd2(i) = hsd1(i)*hsd1p 
!         hsd1p = hsd1(i) 
!         hsd1(i) = hsd1p*(hdi+hdim1) 
!     6   hdim1 = hdi 
!       if (n .eq. 3) hsd1(3) = hsd1(3)+hsd2(2) 
! !                                                                       
! ! test for straight line fit                                            
! !                                                                       
!       con = sumy/float(n) 
!       sum = 0. 
!       do 7 i = 1,n 
!     7   sum = sum+(y(i)-con)**2 
!       if (sum .le. su) go to 23 
! !                                                                       
! ! top of iteration                                                      
! ! cholesky factorization of q*t+h into r                                
! !                                                                       
! !                                                                       
! ! i = 1                                                                 
! !                                                                       
!     8 rd(1) = 1./(q*td(1)+hd(1)) 
!       rnm1(1) = hsd2(1) 
!       yspnm1 = ys(nm1) 
!       rn(1) = q*tsd1(1)+hsd1(1) 
!       yspn = ys(n) 
!       ysp(1) = ys(1) 
!       rsd1i = q*tsd1(2)+hsd1(2) 
!       rsd1(2) = rsd1i*rd(1) 
!       sumnm1 = 0. 
!       sum2 = 0. 
!       sumn = 0. 
!       if (n .eq. 3) go to 11 
!       if (n .eq. 2) go to 12 
! !                                                                       
! ! i = 2                                                                 
! !                                                                       
!       rd(2) = 1./(q*td(2)+hd(2)-rsd1i*rsd1(2)) 
!       rnm1(2) = -rnm1(1)*rsd1(2) 
!       rn(2) = hsd2(2)-rn(1)*rsd1(2) 
!       ysp(2) = ys(2)-rsd1(2)*ysp(1) 
!       if (n .eq. 4) go to 10 
!       do 9 i = 3,nm2 
!         rsd2i = hsd2(i) 
!         rsd1i = q*tsd1(i)+hsd1(i)-rsd2i*rsd1(i-1) 
!         rsd2(i) = rsd2i*rd(i-2) 
!         rsd1(i) = rsd1i*rd(i-1) 
!         rd(i) = 1./(q*td(i)+hd(i)-rsd1i*rsd1(i)                         &
!      &                           -rsd2i*rsd2(i))                        
!         rnm1(i) = -rnm1(i-2)*rsd2(i)-rnm1(i-1)*rsd1(i) 
!         rnm1t = rnm1(i-2)*rd(i-2) 
!         sumnm1 = sumnm1+rnm1t*rnm1(i-2) 
!         rnm1(i-2) = rnm1t 
!         sum2 = sum2+rnm1t*rn(i-2) 
!         yspnm1 = yspnm1-rnm1t*ysp(i-2) 
!         rn(i) = -rn(i-2)*rsd2(i)-rn(i-1)*rsd1(i) 
!         rnt = rn(i-2)*rd(i-2) 
!         sumn = sumn+rnt*rn(i-2) 
!         rn(i-2) = rnt 
!         yspn = yspn-rnt*ysp(i-2) 
!     9   ysp(i) = ys(i)-rsd1(i)*ysp(i-1)-rsd2(i)*ysp(i-2) 
! !                                                                       
! ! i = n-3                                                               
! !                                                                       
!    10 rnm1(nm3) = hsd2(nm1)+rnm1(nm3) 
!       rnm1(nm2) = rnm1(nm2)-hsd2(nm1)*rsd1(nm2) 
!       rnm1t = rnm1(nm3)*rd(nm3) 
!       sumnm1 = sumnm1+rnm1t*rnm1(nm3) 
!       rnm1(nm3) = rnm1t 
!       sum2 = sum2+rnm1t*rn(nm3) 
!       yspnm1 = yspnm1-rnm1t*ysp(nm3) 
!       rnt = rn(nm3)*rd(nm3) 
!       sumn = sumn+rnt*rn(nm3) 
!       rn(nm3) = rnt 
!       yspn = yspn-rnt*ysp(nm3) 
! !                                                                       
! ! i = n-2                                                               
! !                                                                       
!    11 rnm1(nm2) = q*tsd1(nm1)+hsd1(nm1)+rnm1(nm2) 
!       rnm1t = rnm1(nm2)*rd(nm2) 
!       sumnm1 = sumnm1+rnm1t*rnm1(nm2) 
!       rnm1(nm2) = rnm1t 
!       rn(nm2) = hsd2(n)+rn(nm2) 
!       sum2 = sum2+rnm1t*rn(nm2) 
!       yspnm1 = yspnm1-rnm1t*ysp(nm2) 
!       rnt = rn(nm2)*rd(nm2) 
!       sumn = sumn+rnt*rn(nm2) 
!       rn(nm2) = rnt 
!       yspn = yspn-rnt*ysp(nm2) 
! !                                                                       
! ! i = n-1                                                               
! !                                                                       
!    12 rd(nm1) = 1./(q*td(nm1)+hd(nm1)-sumnm1) 
!       ysp(nm1) = yspnm1 
!       rn(nm1) = q*tsd1(n)+hsd1(n)-sum2 
!       rnt = rn(nm1)*rd(nm1) 
!       sumn = sumn+rnt*rn(nm1) 
!       rn(nm1) = rnt 
!       yspn = yspn-rnt*ysp(nm1) 
! !                                                                       
! ! i = n                                                                 
! !                                                                       
!       rdn = q*td(n)+hd(n)-sumn 
!       rd(n) = 0. 
!       if (rdn .gt. 0.) rd(n) = 1./rdn 
!       ysp(n) = yspn 
! !                                                                       
! ! back solve of r(transpose)* r * ysp = ys                              
! !                                                                       
!       ysp(n) = rd(n)*ysp(n) 
!       ysp(nm1) = rd(nm1)*ysp(nm1)-rn(nm1)*ysp(n) 
!       if (n .eq. 2) go to 14 
!       yspn = ysp(n) 
!       yspnm1 = ysp(nm1) 
!       do 13 ibak = 1,nm2 
!         i = nm1-ibak 
!    13   ysp(i) = rd(i)*ysp(i)-rsd1(i+1)*ysp(i+1)                        &
!      &                  -rsd2(i+2)*ysp(i+2)-rnm1(i)*yspnm1              &
!      &                  -rn(i)*yspn                                     
!    14 sum = 0. 
!       delyi1 = (ysp(1)-ysp(n))/(x(1)+p-x(n)) 
!       if (isw .eq. 1) go to 16 
! !                                                                       
! ! calculation of residual norm                                          
! !  - d array                                                            
! !                                                                       
!       do 15 i = 1,nm1 
!         delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i)) 
!         v(i) = (delyi-delyi1)*d(i)*d(i) 
!         sum = sum+v(i)*(delyi-delyi1) 
!    15   delyi1 = delyi 
!       delyi = (ysp(1)-ysp(n))/(x(1)+p-x(n)) 
!       v(n) = (delyi-delyi1)*d(n)*d(n) 
!       go to 18 
! !                                                                       
! ! calculation of residual norm                                          
! !  - d constant                                                         
! !                                                                       
!    16 do 17 i = 1,nm1 
!         delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i)) 
!         v(i) = delyi-delyi1 
!         sum = sum+v(i)*(delyi-delyi1) 
!    17   delyi1 = delyi 
!       delyi = (ysp(1)-ysp(n))/(x(1)+p-x(n)) 
!       v(n) = delyi-delyi1 
!    18   sum = sum+v(n)*(delyi-delyi1) 
! !                                                                       
! ! test for convergence                                                  
! !                                                                       
!       if (sum .le. su .and. sum .ge. sl .and. q .gt. 0.)                &
!      &    go to 21                                                      
! !                                                                       
! ! calculation of newton correction                                      
! !                                                                       
!       f = 0. 
!       g = 0. 
!       rnm1sm = 0. 
!       rnsm = 0. 
!       im1 = n 
!       if (n .eq. 2) go to 20 
!       wim2 = 0. 
!       wim1 = 0. 
!       do 19 i = 1,nm2 
!         tui = tsd1(i)*ysp(im1)+td(i)*ysp(i)                             &
!      &                        +tsd1(i+1)*ysp(i+1)                       
!         wi = tui-rsd1(i)*wim1-rsd2(i)*wim2 
!         rnm1sm = rnm1sm-rnm1(i)*wi 
!         rnsm = rnsm-rn(i)*wi 
!         f = f+tui*ysp(i) 
!         g = g+wi*wi*rd(i) 
!         im1 = i 
!         wim2 = wim1 
!    19   wim1 = wi 
!    20 tui = tsd1(nm1)*ysp(im1)+td(nm1)*ysp(nm1)                         &
!      &      +tsd1(n)*ysp(n)                                             
!       wi = tui+rnm1sm 
!       f = f+tui*ysp(nm1) 
!       g = g+wi*wi*rd(nm1) 
!       tui = tsd1(n)*ysp(nm1)+td(n)*ysp(n)                               &
!      &      +tsd1(1)*ysp(1)                                             
!       wi = tui+rnsm-rn(nm1)*wi 
!       f = f+tui*ysp(n) 
!       g = g+wi*wi*rd(n) 
!       h = f-q*g 
!       if (h .le. 0. .and. q .gt. 0.) go to 21 
! !                                                                       
! ! update q - newton step                                                
! !                                                                       
!       step = (sum-sqrt(sum*sl))/h 
!       if (sl .ne. 0.) step = step*sqrt(sum/sl) 
!       q = q+step 
!       go to 8 
! !                                                                       
! ! store smoothed y-values and second derivatives                        
! !                                                                       
!    21 do 22 i = 1,n 
!         ys(i) = y(i)-v(i) 
!    22   ysp(i) = q*ysp(i) 
!       return 
! !                                                                       
! ! store constant ys and zero ysp                                        
! !                                                                       
!    23 do 24 i = 1,n 
!         ys(i) = con 
!    24   ysp(i) = 0. 
!       return 
! !                                                                       
! ! n less than 2                                                         
! !                                                                       
!    25 ierr = 1 
!       return 
! !                                                                       
! ! s negative                                                            
! !                                                                       
!    26 ierr = 2 
!       return 
! !                                                                       
! ! eps negative or greater than 1                                        
! !                                                                       
!    27 ierr = 3 
!       return 
! !                                                                       
! ! x-values not strictly increasing                                      
! !                                                                       
!    28 ierr = 4 
!       return 
! !                                                                       
! ! weight non-positive                                                   
! !                                                                       
!    29 ierr = 5 
!       return 
! !                                                                       
! ! incorrect period                                                      
! !                                                                       
!    30 ierr = 6 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                         
!       subroutine curvss (n,x,y,d,isw,s,eps,ys,ysp,sigma,td,             &
!      &                   tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,              &
!      &                   ierr)                                          
! !                                                                       
!       integer n,isw,ierr 
!       real x(n),y(n),d(n),s,eps,ys(n),ysp(n),sigma,td(n),               &
!      &     tsd1(n),hd(n),hsd1(n),hsd2(n),rd(n),rsd1(n),                 &
!      &     rsd2(n),v(n)                                                 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine determines the parameters necessary to                
! ! compute a smoothing spline under tension. for a given                 
! ! increasing sequence of abscissae (x(i)), i = 1,..., n and             
! ! associated ordinates (y(i)), i = 1,..., n, the function               
! ! determined minimizes the summation from i = 1 to n-1 of               
! ! the square of the second derivative of f plus sigma                   
! ! squared times the difference of the first derivative of f             
! ! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all               
! ! functions f with two continuous derivatives such that the             
! ! summation of the square of (f(x(i))-y(i))/d(i) is less                
! ! than or equal to a given constant s, where (d(i)), i = 1,             
! ! ..., n are a given set of observation weights. the                    
! ! function determined is a spline under tension with third              
! ! derivative discontinuities at (x(i)), i = 2,..., n-1. for             
! ! actual computation of points on the curve it is necessary             
! ! to call the function curv2.                                           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   n is the number of values to be smoothed (n.ge.2).                  
! !                                                                       
! !   x is an array of the n increasing abscissae of the                  
! !   values to be smoothed.                                              
! !                                                                       
! !   y is an array of the n ordinates of the values to be                
! !   smoothed, (i. e. y(k) is the functional value                       
! !   corresponding to x(k) ).                                            
! !                                                                       
! !   d is a parameter containing the observation weights.                
! !   this may either be an array of length n or a scalar                 
! !   (interpreted as a constant). the value of d                         
! !   corresponding to the observation (x(k),y(k)) should                 
! !   be an approximation to the standard deviation of error.             
! !                                                                       
! !   isw contains a switch indicating whether the parameter              
! !   d is to be considered a vector or a scalar,                         
! !          = 0 if d is an array of length n,                            
! !          = 1 if d is a scalar.                                        
! !                                                                       
! !   s contains the value controlling the smoothing. this                
! !   must be non-negative. for s equal to zero, the                      
! !   subroutine does interpolation, larger values lead to                
! !   smoother funtions. if parameter d contains standard                 
! !   deviation estimates, a reasonable value for s is                    
! !   float(n).                                                           
! !                                                                       
! !   eps contains a tolerance on the relative precision to               
! !   which s is to be interpreted. this must be greater than             
! !   or equal to zero and less than equal or equal to one. a             
! !   reasonable value for eps is sqrt(2./float(n)).                      
! !                                                                       
! !   ys is an array of length at least n.                                
! !                                                                       
! !   ysp is an array of length at least n.                               
! !                                                                       
! !   sigma contains the tension factor. this value indicates             
! !   the degree to which the first derivative part of the                
! !   smoothing functional is emphasized. if sigma is nearly              
! !   zero (e. g. .001) the resulting curve is approximately a            
! !   cubic spline. if sigma is large (e. g. 50.) the                     
! !   resulting curve is nearly a polygonal line. if sigma                
! !   equals zero a cubic spline results. a standard value for            
! !   sigma is approximately 1.                                           
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, and v are                 
! !   arrays of length at least n which are used for scratch              
! !   storage.                                                            
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   ys contains the smoothed ordinate values.                           
! !                                                                       
! !   ysp contains the values of the second derivative of the             
! !   smoothed curve at the given nodes.                                  
! !                                                                       
! !   ierr contains an error flag,                                        
! !        = 0 for normal return,                                         
! !        = 1 if n is less than 2,                                       
! !        = 2 if s is negative,                                          
! !        = 3 if eps is negative or greater than one,                    
! !        = 4 if x-values are not strictly increasing,                   
! !        = 5 if a d-value is non-positive.                              
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n, x, y, d, isw, s, eps, and sigma are unaltered.                   
! !                                                                       
! ! this subroutine references package modules terms and                  
! ! snhcsh.                                                               
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       if (n .lt. 2) go to 16 
!       if (s .lt. 0.) go to 17 
!       if (eps .lt. 0. .or. eps .gt. 1.) go to 18 
!       ierr = 0 
!       p = 0. 
!       v(1) = 0. 
!       v(n) = 0. 
!       ysp(1) = 0. 
!       ysp(n) = 0. 
!       if (n .eq. 2) go to 14 
!       rsd1(1) = 0. 
!       rd(1) = 0. 
!       rsd2(n) = 0. 
!       rdim1 = 0. 
!       yspim2 = 0. 
! !                                                                       
! ! denormalize tension factor                                            
! !                                                                       
!       sigmap = abs(sigma)*float(n-1)/(x(n)-x(1)) 
! !                                                                       
! ! form t matrix and second differences of y into ys                     
! !                                                                       
!       nm1 = n-1 
!       nm3 = n-3 
!       delxi1 = 1. 
!       delyi1 = 0. 
!       dim1 = 0. 
!       do 1 i = 1,nm1 
!         delxi = x(i+1)-x(i) 
!         if (delxi .le. 0.) go to 19 
!         delyi = (y(i+1)-y(i))/delxi 
!         ys(i) = delyi-delyi1 
!         call terms (di,tsd1(i+1),sigmap,delxi) 
!         td(i) = di+dim1 
!         hd(i) = -(1./delxi+1./delxi1) 
!         hsd1(i+1) = 1./delxi 
!         delxi1 = delxi 
!         delyi1 = delyi 
!     1   dim1 = di 
! !                                                                       
! ! calculate lower and upper tolerances                                  
! !                                                                       
!       sl = s*(1.-eps) 
!       su = s*(1.+eps) 
!       if (isw .eq. 1) go to 3 
! !                                                                       
! ! form h matrix - d array                                               
! !                                                                       
!       if (d(1) .le. 0. .or. d(2) .le. 0.) go to 20 
!       betapp = 0. 
!       betap = 0. 
!       alphap = 0. 
!       do 2 i = 2,nm1 
!         alpha = hd(i)*d(i)*d(i) 
!         if (d(i+1) .le. 0.) go to 20 
!         beta = hsd1(i+1)*d(i+1)*d(i+1) 
!         hd(i) = (hsd1(i)*d(i-1))**2+alpha*hd(i)                         &
!      &                             +beta*hsd1(i+1)                      
!         hsd2(i) = hsd1(i)*betapp 
!         hsd1(i) = hsd1(i)*(alpha+alphap) 
!         alphap = alpha 
!         betapp = betap 
!     2   betap = beta 
!       go to 5 
! !                                                                       
! ! form h matrix - d constant                                            
! !                                                                       
!     3 if (d(1) .le. 0.) go to 20 
!       sl = d(1)*d(1)*sl 
!       su = d(1)*d(1)*su 
!       hsd1p = 0. 
!       hdim1 = 0. 
!       do 4 i = 2,nm1 
!         hdi = hd(i) 
!         hd(i) = hsd1(i)*hsd1(i)+hdi*hdi+hsd1(i+1)*hsd1(i+1) 
!         hsd2(i) = hsd1(i)*hsd1p 
!         hsd1p = hsd1(i) 
!         hsd1(i) = hsd1p*(hdi+hdim1) 
!     4   hdim1 = hdi 
! !                                                                       
! ! top of iteration                                                      
! ! cholesky factorization of p*t+h into r                                
! !                                                                       
!     5 do 6 i = 2,nm1 
!         rsd2i = hsd2(i) 
!         rsd1i = p*tsd1(i)+hsd1(i)-rsd2i*rsd1(i-1) 
!         rsd2(i) = rsd2i*rdim1 
!         rdim1 = rd(i-1) 
!         rsd1(i) = rsd1i*rdim1 
!         rd(i) = 1./(p*td(i)+hd(i)-rsd1i*rsd1(i)                         &
!      &                           -rsd2i*rsd2(i))                        
!         ysp(i) = ys(i)-rsd1(i)*ysp(i-1)-rsd2(i)*yspim2 
!     6   yspim2 = ysp(i-1) 
! !                                                                       
! ! back solve of r(transpose)* r * ysp = ys                              
! !                                                                       
!       ysp(nm1) = rd(nm1)*ysp(nm1) 
!       if (n .eq. 3) go to 8 
!       do 7 ibak = 1,nm3 
!         i = nm1-ibak 
!     7   ysp(i) = rd(i)*ysp(i)-rsd1(i+1)*ysp(i+1)                        &
!      &                       -rsd2(i+2)*ysp(i+2)                        
!     8 sum = 0. 
!       delyi1 = 0. 
!       if (isw .eq. 1) go to 10 
! !                                                                       
! ! calculation of residual norm                                          
! !  - d array                                                            
! !                                                                       
!       do 9 i = 1,nm1 
!         delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i)) 
!         v(i) = (delyi-delyi1)*d(i)*d(i) 
!         sum = sum+v(i)*(delyi-delyi1) 
!     9   delyi1 = delyi 
!       v(n) = -delyi1*d(n)*d(n) 
!       go to 12 
! !                                                                       
! ! calculation of residual norm                                          
! !  - d constant                                                         
! !                                                                       
!    10 do 11 i = 1,nm1 
!         delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i)) 
!         v(i) = delyi-delyi1 
!         sum = sum+v(i)*(delyi-delyi1) 
!    11   delyi1 = delyi 
!       v(n) = -delyi1 
!    12 sum = sum-v(n)*delyi1 
! !                                                                       
! ! test for convergence                                                  
! !                                                                       
!       if (sum .le. su) go to 14 
! !                                                                       
! ! calculation of newton correction                                      
! !                                                                       
!       f = 0. 
!       g = 0. 
!       wim2 = 0. 
!       wim1 = 0. 
!       do 13 i = 2,nm1 
!         tui = tsd1(i)*ysp(i-1)+td(i)*ysp(i)                             &
!      &                        +tsd1(i+1)*ysp(i+1)                       
!         wi = tui-rsd1(i)*wim1-rsd2(i)*wim2 
!         f = f+tui*ysp(i) 
!         g = g+wi*wi*rd(i) 
!         wim2 = wim1 
!    13   wim1 = wi 
!       h = f-p*g 
!       if (h .le. 0.) go to 14 
! !                                                                       
! ! update p - newton step                                                
! !                                                                       
!       step = (sum-sqrt(sum*sl))/h 
!       if (sl .ne. 0.) step = step*sqrt(sum/sl) 
!       p = p+step 
!       go to 5 
! !                                                                       
! ! store smoothed y-values and second derivatives                        
! !                                                                       
!    14 do 15 i = 1,n 
!         ys(i) = y(i)-v(i) 
!    15   ysp(i) = p*ysp(i) 
!       return 
! !                                                                       
! ! n less than 2                                                         
! !                                                                       
!    16 ierr = 1 
!       return 
! !                                                                       
! ! s negative                                                            
! !                                                                       
!    17 ierr = 2 
!       return 
! !                                                                       
! ! eps negative or greater than 1                                        
! !                                                                       
!    18 ierr = 3 
!       return 
! !                                                                       
! ! x-values not strictly increasing                                      
! !                                                                       
!    19 ierr = 4 
!       return 
! !                                                                       
! ! weight non-positive                                                   
! !                                                                       
!    20 ierr = 5 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                         
!       function intrvl (t,x,n) 
! !                                                                       
!       integer n 
!       real t,x(n) 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function determines the index of the interval                    
! ! (determined by a given increasing sequence) in which                  
! ! a given value lies.                                                   
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t is the given value.                                               
! !                                                                       
! !   x is a vector of strictly increasing values.                        
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   n is the length of x (n .ge. 2).                                    
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   intrvl returns an integer i such that                               
! !                                                                       
! !          i =  1       if         e   t .lt. x(2)  ,                   
! !          i =  n-1     if x(n-1) .le. t            ,                   
! !          otherwise       x(i)  .le. t .le. x(i+1),                    
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       save i 
!       data i /1/ 
! !                                                                       
!       tt = t 
! !                                                                       
! ! check for illegal i                                                   
! !                                                                       
!       if (i .ge. n) i = n/2 
! !                                                                       
! ! check old interval and extremes                                       
! !                                                                       
!       if (tt .lt. x(i)) then 
!         if (tt .le. x(2)) then 
!           i = 1 
!           intrvl = 1 
!           return 
!         else 
!           il = 2 
!           ih = i 
!         end if 
!       else if (tt .le. x(i+1)) then 
!         intrvl = i 
!         return 
!       else if (tt .ge. x(n-1)) then 
!         i = n-1 
!         intrvl = n-1 
!         return 
!       else 
!         il = i+1 
!         ih = n-1 
!       end if 
! !                                                                       
! ! binary search loop                                                    
! !                                                                       
!     1 i = (il+ih)/2 
!       if (tt .lt. x(i)) then 
!          ih = i 
!       else if (tt .gt. x(i+1)) then 
!          il = i+1 
!       else 
!          intrvl = i 
!          return 
!       end if 
!       go to 1 
!       END function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
!       function intrvp (t,x,n,p,tp) 
! !                                                                       
!       integer n 
!       real t,x(n),p,tp 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this function determines the index of the interval                    
! ! (determined by a given increasing sequence) in which a                
! ! given value lies, after translating the value to within               
! ! the correct period.  it also returns this translated value.           
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   t is the given value.                                               
! !                                                                       
! !   x is a vector of strictly increasing values.                        
! !                                                                       
! !   n is the length of x (n .ge. 2).                                    
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   p contains the period.                                              
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   tp contains a translated value of t (i. e. x(1) .le. tp,            
! !   tp .lt. x(1)+p, and tp = t + k*p for some integer k).               
! !                                                                       
! !   intrvl returns an integer i such that                               
! !                                                                       
! !          i = 1       if             tp .lt. x(2)  ,                   
! !          i = n       if   x(n) .le. tp            ,                   
! !          otherwise       x(i)  .le. tp .lt. x(i+1),                   
! !                                                                       
! ! none of the input parameters are altered.                             
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       save i 
!       data i /1/ 
! !                                                                       
!       nper = (t-x(1))/p 
!       tp = t-float(nper)*p 
!       if (tp .lt. x(1)) tp = tp+p 
!       tt = tp 
! !                                                                       
! ! check for illegal i                                                   
! !                                                                       
!       if (i .ge. n) i = n/2 
! !                                                                       
! ! check old interval and extremes                                       
! !                                                                       
!       if (tt .lt. x(i)) then 
!         if (tt .le. x(2)) then 
!           i = 1 
!           intrvp = 1 
!           return 
!         else 
!           il = 2 
!           ih = i 
!         end if 
!       else if (tt .le. x(i+1)) then 
!         intrvp = i 
!         return 
!       else if (tt .ge. x(n)) then 
!         i = n 
!         intrvp = n 
!         return 
!       else 
!         il = i+1 
!         ih = n 
!       end if 
! !                                                                       
! ! binary search loop                                                    
! !                                                                       
!     1 i = (il+ih)/2 
!       if (tt .lt. x(i)) then 
!          ih = i 
!       else if (tt .gt. x(i+1)) then 
!          il = i+1 
!       else 
!          intrvp = i 
!          return 
!       end if 
!       go to 1 
!       END function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
!                                           
!       subroutine snhcsh2 (sinhm,coshm,x,isw) 
! !                                                                       
!       integer isw 
!       real sinhm,coshm,x 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine returns approximations to                             
! !       sinhm(x) = sinh(x)/x-1                                          
! !       coshm(x) = cosh(x)-1                                            
! ! and                                                                   
! !       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)                             
! ! with relative error less than 1.0e-6                                  
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   x contains the value of the independent variable.                   
! !                                                                       
! !   isw indicates the function desired                                  
! !           = -1 if only sinhm is desired,                              
! !           =  0 if both sinhm and coshm are desired,                   
! !           =  1 if only coshm is desired,                              
! !           =  2 if only coshmm is desired,                             
! !           =  3 if both sinhm and coshmm are desired.                  
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   sinhm contains the value of sinhm(x) if isw .le. 0 or               
! !   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.             
! !   2).                                                                 
! !                                                                       
! !   coshm contains the value of coshm(x) if isw .eq. 0 or               
! !   isw .eq. 1 and contains the value of coshmm(x) if isw               
! !   .ge. 2 (coshm is unaltered if isw .eq. -1).                         
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   x and isw are unaltered.                                            
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       data sp13/.3029390e-5/,                                           &
!      &     sp12/.1975135e-3/,                                           &
!      &     sp11/.8334261e-2/,                                           &
!      &     sp10/.1666665e0/                                             
!       data sp24/.3693467e-7/,                                           &
!      &     sp23/.2459974e-5/,                                           &
!      &     sp22/.2018107e-3/,                                           &
!      &     sp21/.8315072e-2/,                                           &
!      &     sp20/.1667035e0/                                             
!       data sp33/.6666558e-5/,                                           &
!      &     sp32/.6646307e-3/,                                           &
!      &     sp31/.4001477e-1/,                                           &
!      &     sq32/.2037930e-3/,                                           &
!      &     sq31/-.6372739e-1/,                                          &
!      &     sq30/.6017497e1/                                             
!       data sp43/.2311816e-4/,                                           &
!      &     sp42/.2729702e-3/,                                           &
!      &     sp41/.9868757e-1/,                                           &
!      &     sq42/.1776637e-3/,                                           &
!      &     sq41/-.7549779e-1/,                                          &
!      &     sq40/.9110034e1/                                             
!       data cp4/.2982628e-6/,                                            &
!      &     cp3/.2472673e-4/,                                            &
!      &     cp2/.1388967e-2/,                                            &
!      &     cp1/.4166665e-1/,                                            &
!      &     cp0/.5000000e0/                                              
! !                                                                       
!       ax = abs(x) 
!       if (isw .ge. 0) go to 5 
! !                                                                       
! ! sinhm approximation                                                   
! !                                                                       
!       if (ax .gt. 4.45) go to 2 
!       xs = ax*ax 
!       if (ax .gt. 2.3) go to 1 
! !                                                                       
! ! sinhm approximation on (0.,2.3)                                       
! !                                                                       
!       sinhm = xs*(((sp13*xs+sp12)*xs+sp11)*xs+sp10) 
!       return 
! !                                                                       
! ! sinhm approximation on (2.3,4.45)                                     
! !                                                                       
!     1 sinhm = xs*((((sp24*xs+sp23)*xs+sp22)*xs+sp21)                    &
!      &               *xs+sp20)                                          
!       return 
!     2 if (ax .gt. 7.65) go to 3 
! !                                                                       
! ! sinhm approximation on (4.45,7.65)                                    
! !                                                                       
!       xs = ax*ax 
!       sinhm = xs*(((sp33*xs+sp32)*xs+sp31)*xs+1.)/                      &
!      &             ((sq32*xs+sq31)*xs+sq30)                             
!       return 
!     3 if (ax .gt. 10.1) go to 4 
! !                                                                       
! ! sinhm approximation on (7.65,10.1)                                    
! !                                                                       
!       xs = ax*ax 
!       sinhm = xs*(((sp43*xs+sp42)*xs+sp41)*xs+1.)/                      &
!      &             ((sq42*xs+sq41)*xs+sq40)                             
!       return 
! !                                                                       
! ! sinhm approximation above 10.1                                        
! !                                                                       
!     4 sinhm = exp(ax)/(ax+ax)-1. 
!       return 
! !                                                                       
! ! coshm and (possibly) sinhm approximation                              
! !                                                                       
!     5 if (isw .ge. 2) go to 7 
!       if (ax .gt. 2.3) go to 6 
!       xs = ax*ax 
!       coshm = xs*((((cp4*xs+cp3)*xs+cp2)*xs+cp1)*xs+cp0) 
!       if (isw .eq. 0) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)              &
!      &                              *xs+sp10)                           
!       return 
!     6 expx = exp(ax) 
!       coshm = (expx+1./expx)/2.-1. 
!       if (isw .eq. 0) sinhm = (expx-1./expx)/(ax+ax)-1. 
!       return 
! !                                                                       
! ! coshmm and (possibly) sinhm approximation                             
! !                                                                       
!     7 xs = ax*ax 
!       if (ax .gt. 2.3) go to 8 
!       coshm = xs*(((cp4*xs+cp3)*xs+cp2)*xs+cp1) 
!       if (isw .eq. 3) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)              &
!      &                              *xs+sp10)                           
!       return 
!     8 expx = exp(ax) 
!       coshm = ((expx+1./expx-xs)/2.-1.)/xs 
!       if (isw .eq. 3) sinhm = (expx-1./expx)/(ax+ax)-1. 
!       return 
!       END subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
! 
! 
! 
! 
!       subroutine snhcsh (sinhm,coshm,x,isw) 
! !                                                                       
!       integer isw 
!       real sinhm,coshm,x 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine returns approximations to                             
! !       sinhm(x) = sinh(x)/x-1                                          
! !       coshm(x) = cosh(x)-1                                            
! ! and                                                                   
! !       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)                             
! ! with relative error less than 4.0e-14.                                
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   x contains the value of the independent variable.                   
! !                                                                       
! !   isw indicates the function desired                                  
! !           = -1 if only sinhm is desired,                              
! !           =  0 if both sinhm and coshm are desired,                   
! !           =  1 if only coshm is desired,                              
! !           =  2 if only coshmm is desired,                             
! !           =  3 if both sinhm and coshmm are desired.                  
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !   sinhm contains the value of sinhm(x) if isw .le. 0 or               
! !   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.             
! !   2).                                                                 
! !                                                                       
! !   coshm contains the value of coshm(x) if isw .eq. 0 or               
! !   isw .eq. 1 and contains the value of coshmm(x) if isw               
! !   .ge. 2 (coshm is unaltered if isw .eq. -1).                         
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   x and isw are unaltered.                                            
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       data sp14/.227581660976348e-7/,                                   &
!      &     sp13/.612189863171694e-5/,                                   &
!      &     sp12/.715314759211209e-3/,                                   &
!      &     sp11/.398088289992973e-1/,                                   &
!      &     sq12/.206382701413725e-3/,                                   &
!      &     sq11/-.611470260009508e-1/,                                  &
!      &     sq10/.599999999999986e+1/                                    
!       data sp25/.129094158037272e-9/,                                   &
!      &     sp24/.473731823101666e-7/,                                   &
!      &     sp23/.849213455598455e-5/,                                   &
!      &     sp22/.833264803327242e-3/,                                   &
!      &     sp21/.425024142813226e-1/,                                   &
!      &     sq22/.106008515744821e-3/,                                   &
!      &     sq21/-.449855169512505e-1/,                                  &
!      &     sq20/.600000000268619e+1/                                    
!       data sp35/.155193945864942e-9/,                                   &
!      &     sp34/.511529451668737e-7/,                                   &
!      &     sp33/.884775635776784e-5/,                                   &
!      &     sp32/.850447617691392e-3/,                                   &
!      &     sp31/.428888148791777e-1/,                                   &
!      &     sq32/.933128831061610e-4/,                                   &
!      &     sq31/-.426677570538507e-1/,                                  &
!      &     sq30/.600000145086489e+1/                                    
!       data sp45/.188070632058331e-9/,                                   &
!      &     sp44/.545792817714192e-7/,                                   &
!      &     sp43/.920119535795222e-5/,                                   &
!      &     sp42/.866559391672985e-3/,                                   &
!      &     sp41/.432535234960858e-1/,                                   &
!      &     sq42/.824891748820670e-4/,                                   &
!      &     sq41/-.404938841672262e-1/,                                  &
!      &     sq40/.600005006283834e+1/                                    
!       data cp5/.552200614584744e-9/,                                    &
!      &     cp4/.181666923620944e-6/,                                    &
!      &     cp3/.270540125846525e-4/,                                    &
!      &     cp2/.206270719503934e-2/,                                    &
!      &     cp1/.744437205569040e-1/,                                    &
!      &     cq2/.514609638642689e-4/,                                    &
!      &     cq1/-.177792255528382e-1/,                                   &
!      &     cq0/.200000000000000e+1/                                     
!       data zp4/.664418805876835e-8/,                                    &
!      &     zp3/.218274535686385e-5/,                                    &
!      &     zp2/.324851059327161e-3/,                                    &
!      &     zp1/.244515150174258e-1/,                                    &
!      &     zq2/.616165782306621e-3/,                                    &
!      &     zq1/-.213163639579425e0/,                                    &
!      &     zq0/.240000000000000e+2/                                     
! !                                                                       
!       ax = abs(x) 
!       if (isw .ge. 0) go to 5 
! !                                                                       
! ! sinhm approximation                                                   
! !                                                                       
!       if (ax .gt. 3.9) go to 2 
!       xs = ax*ax 
!       if (ax .gt. 2.2) go to 1 
! !                                                                       
! ! sinhm approximation on (0.,2.2)                                       
! !                                                                       
!       sinhm = xs*((((sp14*xs+sp13)*xs+sp12)*xs+sp11)*xs+1.)/            &
!      &             ((sq12*xs+sq11)*xs+sq10)                             
!       return 
! !                                                                       
! ! sinhm approximation on (2.2,3.9)                                      
! !                                                                       
!     1 sinhm = xs*(((((sp25*xs+sp24)*xs+sp23)*xs+sp22)*xs+sp21)          &
!      &        *xs+1.)/((sq22*xs+sq21)*xs+sq20)                          
!       return 
!     2 if (ax .gt. 5.1) go to 3 
! !                                                                       
! ! sinhm approximation on (3.9,5.1)                                      
! !                                                                       
!       xs = ax*ax 
!       sinhm = xs*(((((sp35*xs+sp34)*xs+sp33)*xs+sp32)*xs+sp31)          &
!      &        *xs+1.)/((sq32*xs+sq31)*xs+sq30)                          
!       return 
!     3 if (ax .gt. 6.1) go to 4 
! !                                                                       
! ! sinhm approximation on (5.1,6.1)                                      
! !                                                                       
!       xs = ax*ax 
!       sinhm = xs*(((((sp45*xs+sp44)*xs+sp43)*xs+sp42)*xs+sp41)          &
!      &        *xs+1.)/((sq42*xs+sq41)*xs+sq40)                          
!       return 
! !                                                                       
! ! sinhm approximation above 6.1                                         
! !                                                                       
!     4 expx = exp(ax) 
!       sinhm = (expx-1./expx)/(ax+ax)-1. 
!       return 
! !                                                                       
! ! coshm and (possibly) sinhm approximation                              
! !                                                                       
!     5 if (isw .ge. 2) go to 7 
!       if (ax .gt. 2.2) go to 6 
!       xs = ax*ax 
!       coshm = xs*(((((cp5*xs+cp4)*xs+cp3)*xs+cp2)*xs+cp1)               &
!      &        *xs+1.)/((cq2*xs+cq1)*xs+cq0)                             
!       if (isw .eq. 0) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)             &
!      &          *xs+sp11)*xs+1.)/((sq12*xs+sq11)*xs+sq10)               
!       return 
!     6 expx = exp(ax) 
!       coshm = (expx+1./expx)/2.-1. 
!       if (isw .eq. 0) sinhm = (expx-1./expx)/(ax+ax)-1. 
!       return 
! !                                                                       
! ! coshmm and (possibly) sinhm approximation                             
! !                                                                       
!     7 xs = ax*ax 
!       if (ax .gt. 2.2) go to 8 
!       coshm = xs*((((zp4*xs+zp3)*xs+zp2)*xs+zp1)*xs+1.)/                &
!      &             ((zq2*xs+zq1)*xs+zq0)                                
!       if (isw .eq. 3) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)             &
!      &          *xs+sp11)*xs+1.)/((sq12*xs+sq11)*xs+sq10)               
!       return 
!     8 expx = exp(ax) 
!       coshm = ((expx+1./expx-xs)/2.-1.)/xs 
!       if (isw .eq. 3) sinhm = (expx-1./expx)/(ax+ax)-1. 
!       return 
!       END SUBROUTINE   
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! 
!                                         
!       subroutine terms (diag,sdiag,sigma,del) 
! !                                                                       
!       real diag,sdiag,sigma,del 
! !                                                                       
! !                                 coded by alan kaylor cline            
! !                           from fitpack -- january 26, 1987            
! !                        a curve and surface fitting package            
! !                      a product of pleasant valley software            
! !                  8603 altus cove, austin, texas 78759, usa            
! !                                                                       
! ! this subroutine computes the diagonal and superdiagonal               
! ! terms of the tridiagonal linear system associated with                
! ! spline under tension interpolation.                                   
! !                                                                       
! ! on input--                                                            
! !                                                                       
! !   sigma contains the tension factor.                                  
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   del contains the step size.                                         
! !                                                                       
! ! on output--                                                           
! !                                                                       
! !                sigma*del*cosh(sigma*del) - sinh(sigma*del)            
! !   diag = del*--------------------------------------------.            
! !                     (sigma*del)**2 * sinh(sigma*del)                  
! !                                                                       
! !                   sinh(sigma*del) - sigma*del                         
! !   sdiag = del*----------------------------------.                     
! !                (sigma*del)**2 * sinh(sigma*del)                       
! !                                                                       
! ! and                                                                   
! !                                                                       
! !   sigma and del are unaltered.                                        
! !                                                                       
! ! this subroutine references package module snhcsh.                     
! !                                                                       
! !-----------------------------------------------------------            
! !                                                                       
!       if (sigma .ne. 0.) go to 1 
!       diag = del/3. 
!       sdiag = del/6. 
!       return 
!     1 sigdel = sigma*del 
!       call snhcsh (sinhm,coshm,sigdel,0) 
!       denom = sigma*sigdel*(1.+sinhm) 
!       diag = (coshm-sinhm)/denom 
!       sdiag = sinhm/denom 
!       return 
!       END subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!



end module                                           