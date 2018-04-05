"Rad.cov" <-
function (x1, x2, p = 1, with.log = TRUE, with.constant = TRUE, 
            C = NA,marginal=FALSE) 
{
#
# mth order thin plate spline radial basis functions
# in d dimensions
# usually called with p 2m-d

#  marginal dummy argument 
#  this should only be called within predict.se.Krig
#  and provides the correct calculation. Because this is 
#  a generalized covariance the marginal variance is not really 
#  defined. 
#

    if( marginal){
    return( rep( 0, nrow( x1)) )}


        if (!is.matrix(x1)) 
            x1 <- as.matrix(x1)
        if (!is.matrix(x2)) 
            x2 <- as.matrix(x2)
        d <- ncol(x1)
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        m <- (d + p)/2
        par <- c(p/2, ifelse((d%%2 == 0) & (with.log), 1, 0))

# compute matrix in FORTRAN
   
 if (is.na(C[1])) {
        temp <- .Fortran("radbas", nd = as.integer(d), x1 = as.double(x1), 
            n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
            par = as.double(par), k = as.double(rep(0, n1 * n2)), 
            PACKAGE="fields")

        temp<- matrix(temp$k, ncol = n2, nrow = n1)
    }

#
# do cross covariance matrix multiplication in FORTRAN
    else {

        temp<- .Fortran("multrb", nd = as.integer(d), x1 = as.double(x1),
              n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2),
              par = as.double(par), 
              c = as.double(C), h = as.double(rep(0, n1)), 
              work = as.double(rep(0, n2)), PACKAGE="fields")$h

    }

#
# multiply by constant if requested

        if (with.constant) {
            temp <- radbas.constant(m, d)*temp
        }

   return( temp)

}
