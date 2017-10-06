Exp.simple.cov<- function (x1, x2, theta = 1, C=NA, marginal=FALSE) 
{

# different locations are the different rows of x1 and x2. 

# this function can return three different results 
# depending on the values of C and marginal. 
# The three cases:
# 1) cross covaraince matrix
# 2) cross covariance matrix times a vector (C)
# 3) the diagonal elements of covariance matrix at locations x1.

# CASE 1:
 if( is.na( C[1])& !marginal){
# rdist finds the cross distance matrix between the
# locations at x1, x2. 
#
  return(  exp(-rdist(x1, x2)/theta) )
  }
  
# CASE 2:
# or return  multiplication of cov( x2,x1) with vector C
  if(!is.na(C[1])){
   return(  exp(-rdist(x1, x2)/theta)%*%C  )
#
# if the rows of X1 are large 
# this line could be replaced by a call to C or FORTRAN
# to make the multiply use less memory.
#
# there are also other algorithms for fast multiplies when 
# X2 is on a grid.  
#
  }

#  CASE 3
# return marginal variance (in this case it is trivial a constant vector 
# with 1.0)
  if( marginal){
   return( rep( 1, nrow( x1)))
  }

}

