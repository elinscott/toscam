"stationary.cov" <-
 function (x1, x2, Covariance="Exponential", Distance="rdist",
               Dist.args=NULL, theta=1.0,C=NA, marginal=FALSE,...)
{

# get covariance function arguments from call

 Cov.args<- list(...)

# coerce x1 and x2 to matrices

    if( is.data.frame( x1)) x1<- as.matrix(x1)

    if (!is.matrix(x1)) 
        x1 <- matrix(c(x1), ncol=1)

    if (missing(x2)) 
        x2 <- x1

    if( is.data.frame( x2)) x2<- as.matrix(x1)

    if (!is.matrix(x2)) 
        x2 <- matrix(c(x2), ncol=1)

    if (length(theta) == 1) 
        theta <- rep(theta, ncol(x1))

# handle special case of 1-d
    if( ncol(x1)==1) { theta<- matrix( c(theta),1,1)}

# handle special case of just diagonal elements of  theta
    if (is.vector(theta)) 
        theta <- diag(theta)
#
# following now treats theta as a full matrix for scaling and rotation. 
#
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    x1 <- x1 %*% t(solve(theta))
    x2 <- x2 %*% t(solve(theta))
#
# locations are now scaled and rotated correctly 
# now apply covariance function to pairwise distance matrix, or multiply 
# by C vector or just find marginal variance

# this if block finds the cross covariance matrix
 if( is.na( C[1]) & !marginal ){
    bigD<- do.call(Distance,c(list(x1=x1,x2=x2), Dist.args))
    return( do.call(Covariance,  c(list(d=bigD),Cov.args)) )}


# or multiply cross covariance by C
 if( !is.na( C[1]) ){
#
# as coded below this is not particularly efficient of memory
# 
   bigD<- do.call(Distance, c(list(x1=x1,x2=x2), Dist.args) )
   return( 
          do.call(Covariance, c( list(d=bigD),Cov.args) ) %*%C )
 }


# or find marginal variance and return  a vector. 
  if( marginal){
    sigma2<- do.call(Covariance, c(list(d=0),Cov.args) )
    return( rep( sigma2, nrow( x1)))}  

# should not get here

}

