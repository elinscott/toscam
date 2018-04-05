"fields.diagonalize" <-
function(A,B){

eigen( A,symmetric=TRUE)-> hold

# square root of A

hold2 <- (t(hold$vectors) * sqrt( 1/hold$values))

#
# note: interated transposes are a quick way to post multiply by a
# diagonal matrix 
#
#
# A.inv.sqrt = hold2
# A.inv = hold%*% t(hold2) 
#
# eigen decomp of  A.inv.sqrt B t( A.inv.sqrt)
# 
 
eigen(  (hold2)%*%B%*%t(hold2) ,symmetric=TRUE)-> hold

# the magic G matrix used throughout fields. 

G<-  t(hold2)%*%hold$vectors
#
# Note:
# G simultaneously diagonalizes two matrices:
# 
# G^T A G= I
# G^T B G= D
#
# and in terms of application we also have the useful 
# diagonalization
#
#  (A +lambda B)^{-1} =  G( I + lambda D)^{-1} G^T

list( G=G, D=hold$values)
}

