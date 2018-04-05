

      subroutine mkpoly(m,dim,npoly)
      integer m,dim
c
c  Purpose: compute the binomial coefficient of m + dim - 1 choose dim.
c  	This is the dimension of the space of polynomials which are in
c  	the null space of the smoothing penalty. Uses Andy Jaworski's
c	binomial coefficient algorithm that only requires integer
c	arithmetic.
c
c  On Entry:
c   m			order of derivatives in the penalty
c   dim	 		dimension of the variables to be splined
c
c  On Exit:
c   mkploy		(m + dim - 1) choose dim
c
c $Header: /fs/image/home/thoar/CVS.REPOS/fields/src/mkpoly.f,v 1.1.1.1 2002/12/04 22:46:16 thoar Exp $
c
      integer i,j,k,k1,kcoef,n, npoly
c 			compute binomial coefficient
c			m + dim - 1 choose dim
      n = m + dim - 1
      k1 = dim
      if (k1 .gt. n .or. k1 .lt. 0) then
         npoly = 0
         return
      endif
      k = k1
      if ((n - k1) .lt. k) then
         k = n - k1
      endif
      kcoef = 1
      j = n - k
      do 10 i = 1, k
         j = j + 1
         kcoef = (kcoef * j) / i
   10 continue
      npoly = kcoef
      return
      end
