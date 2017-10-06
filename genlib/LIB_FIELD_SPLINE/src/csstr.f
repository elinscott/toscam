  



** finds value of h minimizing the  generalized cross-validation
      subroutine csstr(h,nobs,x,y,wght,c,offset,trace,vlam,work,ierr)

      parameter(mxM=20000)
      implicit double precision (a-h,o-z)
         
      REAL*8 h,trace, vlam,c,offset
      REAL*8 x(nobs),y(nobs),wght(nobs)
      REAL*8 work(nobs),diag(mxM),dumm1(1),dumm2(1)
      integer job(3),ideriv,ierr

      data ideriv/0/
       job(1)=3
       job(2)=0
       job(3)=0
       diag(1)=c
       diag(2)=offset
      call css(h,nobs,x,y,wght,work,trace,diag,vlam,ndum,dumm1,dumm2,
     -            job,ideriv,ierr)

      return
      end
