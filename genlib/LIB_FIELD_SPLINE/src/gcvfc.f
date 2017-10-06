         



** finds value of h minimizing the  generalized cross-validation

      double precision function gcvfc(h,nobs,x,y,wght,c,offset,trace)

      parameter(mxM=20000)
      implicit double precision (a-h,o-z)
         
      REAL*8 h,trace
      REAL*8 x(nobs),y(nobs),wght(nobs)
      REAL*8 sy(mxM),diag(mxM),dumm1(1),dumm2(1)
      integer job(3),ideriv

      data ideriv/0/
      data rinf/1e20/
       job(1)=3
       job(2)=0
       job(3)=0
       diag(1)=c
       diag(2)=offset
      call css(h,nobs,x,y,wght,sy,trace,diag,vlam,ndum,dumm1,dumm2,
     -            job,ideriv,ierr)
      gcvfc= vlam
c
c      write(*, 100) h, c, offset, trace, vlam
 100   format( 5e12.4)
c
c
      return
      end
