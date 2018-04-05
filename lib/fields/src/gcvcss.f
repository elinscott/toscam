
      subroutine gcvcss(n,x,y,wt,c,offset,nstep,maxit,
     - hmin, hmax,hopt,vopt,
     -     tropt,mxstep,fout,
     -                 ierr)
      implicit double precision (a-h,o-z)
      real*8  fout(mxstep,3)
      real*8  x(n),y(n),wt(n)
      real*8  hopt,vopt,trace, hmin,hmax
      
      data tau,tausq/.6180339d0,.3819660d0/
**** coarse search in bandwidth h: this search feeds (hl,hr) to golden search
c compute range for coarse search   
      if (mxstep .lt. nstep) then
         ierr=10
         return
      endif
c
C  rough idea of lambda serach grid
c
      if( hmin.ge.hmax) then 
      xmin= x(1)
      xmax= x(1)
      do 500 i=1,n
          if( x(i).lt.xmin) xmin=x(i)
          if( x(i).gt.xmax) xmax=x(i)
 500  continue
      range = xmax- xmin

      xnobs=n
c**  c= cost value on GCV function there is a pole at trp=2+ nobs/c
      trp=(Xnobs-offset)*.97/c 
      if( trp.le.2.0)  then
c** offset is too big!
            ierr=11
       return
       endif
      hmin=   -0.755762d0 +  0.706693d0*log(xnobs)
     *  +  0.01884722d0*(log(xnobs)**2)  -4.918114d0*log(trp)
     *  +  0.0879931d0*(log(trp)**2)
c      write(*,*) hmin, range, xnobs 
      hmin= hmin + log(range)*3.d0-log(xnobs)

      trp=2.01d0
      hmax=   -0.755762d0 +  0.706693d0*log(xnobs)
     *  +  0.01884722d0*(log(xnobs)**2)  -4.918114d0*log(trp)
     *  +  0.0879931d0*(log(trp)**2)
 
      hmax= hmax + log(range)*3.d0-log(xnobs)
      endif

      hstep=(hmax-hmin)/(nstep-1)
      do 3 j=0,nstep-1
         h=hmin+j*hstep
         gcvh= gcvfc(h,n,x,y,wt,c,offset,trace)
         fout(j+1,1)=h
         fout(j+1,2)=trace
         fout(j+1,3)=gcvh
         if ((gcvh .lt. gcvmin) .or. (j .eq. 0)) then
            hopt=h
            best=h
            gcvmin=gcvh
            trbest=trace
         endif
c          write(*,5000) h, trace,gcvh
 5000   format(5e12.4)
c
  3   continue
c
c      write(*,*) ' Crude search of optimal h completed'
c
**** fast return if crude search minimum gcv is hmin or hmax
           hopt= best
           vopt=gcvmin
           tropt=trbest
      if( (best.eq.hmin) .or. (best.eq.h) ) then
           ierr=-1
           return
      endif


c**** start values for golden search
      hl=best-hstep
      hr=best+hstep

*** Golden section search for min gcv on interval (hl,hr)--maxit(=15) iterations
***   gcv(h) must be quasiconvex on initial interval (hl,hr).
***   On return, h=abscissa of minimum, v=gcv(h)=actual minimum achieved.
***   Interval of uncertainty is (hl,hr), hl <= hlm <= hrm <= hr.

      gcvhl=gcvfc(hl,n,x,y,wt,c,offset,trace)
      gcvhr=gcvfc(hr,n,x,y,wt,c,offset,trace)
      hlm=hl*tau+hr*tausq
      hrm=hl+hr-hlm
      gcvhlm=gcvfc(hlm,n,x,y,wt,c,offset,trchlm)
      gcvhrm=gcvfc(hrm,n,x,y,wt,c,offset,trchrm)

      do 5 it=1,maxit
         if( gcvhlm .ge. gcvhrm ) then
            if (gcvhl .lt. gcvhlm) then
               err= gcvhl/gcvhlm
c              write(*,10) it,err,hl,hlm
               ierr=2
               return
            endif
c
            hl=hlm
            gcvhl=gcvhlm
            hlm=hrm
            hrm=hrm+(hrm-hl)*tau
            gcvhlm=gcvhrm
            gcvhrm=gcvfc(hrm,n,x,y,wt,c,offset,trchrm) 
         else
            if (gcvhr .lt. gcvhrm) then
               err= gcvhrm/gcvhr
c              write(*,11) it,err,hr,hrm
               ierr=2
               return
            endif
c
            hr=hrm
            gcvhr=gcvhrm
            hrm=hlm
            hlm=hlm+(hlm-hr)*tau
            gcvhrm=gcvhlm
            gcvhlm=gcvfc(hlm,n,x,y, wt,c,offset,trchlm)
         endif
5     continue

**** finished -- take the best h so far
      if( gcvhlm .ge. gcvhrm) then
         hopt=hrm
         vopt=gcvhrm
         tropt=trchrm
      else
         hopt=hlm
         vopt=gcvhlm
         tropt=trchlm
      endif
c
10    format(' gcv(h) is NOT quasiconvex ',/,
     - ' iter,error,hl,hlm: ',I4,e15.5,2e15.5)
11    format(' gcv(h) is NOT quasiconvex ',/,
     - ' iter,error,hr,hrm: ',I4,e15.5,2e15.5)
c
c
c     write(*,*) 'hopt', hopt
      return
      end
