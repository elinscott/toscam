"Exp.cov" <-
function (x1, x2, theta = rep(1, ncol(x1)), p = 1, C = NA,marginal=FALSE) 
{
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    if (length(theta) == 1) 
        theta <- rep(theta, ncol(x1))
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
# scale the coordinates by theta
# a more general scaling by a matrix is done in stationary.cov
    x1 <- scale(x1 , center=FALSE, scale=theta)
    x2 <- scale(x2 , center=FALSE, scale=theta)
    par <- p
#
# there are three possible actions listed below:

# find cross covariance matrix
    if (is.na(C[1])& !marginal) {
       return( exp(-rdist(x1, x2)^p) )
    }
#
# multiply cross covariance matrix by C
# in this case implemented in FORTRAN
#
    if(!is.na( C[1])) {
       return(
        .Fortran("multeb", nd = as.integer(d), x1 = as.double(x1), 
            n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
            par = as.double(par), c = as.double(C), h = as.double(rep(0, 
                n1)), work = as.double(rep(0, n2)),PACKAGE="fields")$h
      )
    }
#
# return marginal variance ( 1.0 in this case)

    if( marginal){
    return( rep( 1, nrow(x1)) )}

}
