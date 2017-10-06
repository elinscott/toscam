"splint" <-
function (x, y, xgrid, derivative = 0) 
{
    if (is.matrix(x)) {
      if( ncol(x)>1){
        xgrid <- y
        y <- x[, 2]
        x <- x[, 1]}
    }
    if (is.list(x)) {
        xgrid <- y
        y <- x$y
        x <- x$x
    }
    ind <- !duplicated(x)
    x <- x[ind]
    y <- y[ind]
    if ((derivative > 2) | (derivative < 0)) 
        stop("derivative must be 0,1,2")
    if (length(x) != length(y)) 
        stop("Lengths of x and y must match")
    n <- length(x)
    .Fortran("css", as.double(0), as.integer(n), as.double(x), 
        as.double(y), as.double(rep(0, n)), as.double(rep(1, 
            n)), as.double(1), as.double(1), as.double(1), as.integer(length(xgrid)), 
        as.double(xgrid), ygrid = as.double(rep(0, length(xgrid))), 
        as.integer(c(2, 3, 0)), as.integer(derivative), 
        as.integer(0),PACKAGE="fields")$ygrid
}
