"rdist" <-
function (x1, x2) 
{
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if( missing( x2)) x2<- x1
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    par <- c(1/2, 0)
    matrix((.Fortran("radbas", nd = as.integer(d), x1 = as.double(x1), 
        n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
        par = as.double(par), k = as.double(rep(0, n1 * n2)),PACKAGE="fields")$k), 
        ncol = n2, nrow = n1)
}
