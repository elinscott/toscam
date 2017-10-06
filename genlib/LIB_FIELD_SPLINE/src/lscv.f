
      subroutine lscv(nh,h,n,m,x,cvf)
c*** computes  computes a least squares cross validation function 
C**** for a density estimate using a normal kernel
c**** 
      REAL*8 x(n,m),h(nh),cvf(nh),dist
      real*8 hm,sum
      real*8 rh, con1,con2,k10,k20,ks0
      real*8 tk,tk2
      integer n,nh,m, j,jj,jm1,k,ih

      con1= .39894228
      k10= con1**m
      con2= .39894228/sqrt(2.0)
      k20= con2**m
      ks0=   k20- 2.0* k10


      do 100 ih=1, nh
       rh=1/h(ih)**2
      hm= h(ih)**m
      sum =0
       
c**** compute off diagonal elements of integral of fhat squared
      DO 2 J=2,N              
            jm1= j-1
       do 3 k=1,jm1
          dist=0.0
          do 6 jj=1,m
              dist= dist +  ( x(j,jj)-x(k,jj))**2
  6       continue
          tk2=  exp( -.25*dist*rh)
          tk= tk2*tk2
         
         sum=sum + k20*tk2 - 2*k10*tk
    3 continue

    2 continue
c**** add together off diag and diagonal terms
c**** this is proportional to the integral of fhat **2
       sum = 2*sum + n* ks0
c**** add in contribution from cross term   integral fhat*f
       sum = ( sum / (n*n) + 2*k10/(n) )/hm
       cvf(ih)=sum
  100 continue


      return    
      END    
