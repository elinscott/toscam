

       subroutine m2deb( nd,x1,n1, x2,n2, par, c,nc,h,work)
       implicit double precision (a-h,o-z)
       integer n1,n2, ir,nc,nd,ic
       
       real*8 par(nd),x1(n1,nd), x2(n2,nd), c(n2,nc), h(n1,nc)
       real*8 work( nc),  xtemp1, xtemp2,temp, temp2

c****** work aray must be dimensioned to size n2
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce paging 
       do 4 ic=1,nc
            work(ic)=0
 4     continue
       do 5 ir= 1, n1
       xtemp1= x1(ir,1)
       xtemp2= x1(ir,2)

c evaluate row of exp covaraince matrix       
       
       do 10 j =1,n2

        temp= sqrt((xtemp1- x2(j,1))**2 +(xtemp2- x2(j,2))**2)
        temp2= exp( -temp)
        do 15 ic= 1,nc
        work(ic)= work(ic)+  temp2 * c(j,ic)
 15     continue

 10    continue
        
          do 25 ic=1,nc
          h(ir,ic)= work(ic)
          work( ic)=0
 25        continue
 5      continue

       return
       end
