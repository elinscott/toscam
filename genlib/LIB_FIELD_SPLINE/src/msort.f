
         subroutine msort( m,n,x)
c sorts the columns of mXn matrix in place
         real*8 x(m,1)
         integer m,n,j
         do 10 j=1,n
                call hsort(x(1,j),m)
 10      continue
         return
         end
