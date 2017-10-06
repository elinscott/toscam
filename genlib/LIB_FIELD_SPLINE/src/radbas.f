c**** subroutine to fill in the omega ( or K) matrix for 
c**** ridge regression S funcion
c**** K_ij= radfun( distance( x1_i, x2_j))
c
       subroutine radbas( nd,x1,n1, x2,n2, par, k)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,ic
       
       real*8 par(1),x1(n1,nd), x2(n2,nd), k(n1,n2)
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce paging 

       do 5 ic= 1, nd
       do 10 j =1,n2
            xtemp= x2(j,ic)
            do  15 i= 1, n1
c
c** accumulate sqared differences
c 
               k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
 15             continue
 10    continue
 5      continue
c**** at this point k( i,j) is the squared distnace between x1_i and x2_j
c*** now evaluate radial basis functions
         nbig= n1*n2
c***** Now evalute the radial basis functions with the
c      distances. radfun will just loop through the matrix 
c as stacked column vectors. 

         call radfun( nbig,k(1,1),par)

       return
       end
