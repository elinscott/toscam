




c**** subroutine to fill in the omega ( or K) matrix for 
c**** GASP type covariance
c**** K_ij= exp(-( sum_k( |x1_ik- x2_jk|**par(k)) )
C**** 
C*** in the code k is really ic
c
       subroutine gaspbs( nd,x1,n1, x2,n2, par, k)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,ic
       
       real*8 par(nd),x1(n1,nd), x2(n2,nd), k(n1,n2)
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce paging 

       do 5 ic= 1, nd
       do 10 j =1,n2
            xtemp= x2(j,ic)
            do  15 i= 1, n1
c
c** pth power differences

               k(i,j)=  (abs(x1(i,ic)- xtemp))**par(ic) + k(i,j)
 15             continue
 10    continue
 5      continue
c**** at this point k( i,j) is the 
c***** sum of differences each raised to the par(ic) power
c***
c*** now evaluate like radial basis functions
         nbig= n1*n2
c***** Now evalute the radial basis functions with the
c      distances. gaspf will just loop through the matrix 
c as stacked column vectors. 

         call gaspfn( nbig,k(1,1),par)

       return
       end
