"fields.mkpoly" <-
function (x, m = 2) 
{

 if (m < 1) stop("'m' has to be larger than zero.")
 if (!is.matrix(x)) x <- as.matrix(x)	    

    d <- ncol(x)
    n <- nrow(x)
    nterms <- .Fortran("mkpoly", as.integer(m), as.integer(d), 
        nterms = as.integer(0), PACKAGE = "fields")$nterms
    temp <- .Fortran("dmaket", m = as.integer(m), n = as.integer(n), 
        dim = as.integer(d), des = as.double(x), lddes = as.integer(n), 
        npoly = as.integer(nterms), tmatrix = as.double(rep(0, 
            n * (nterms))), ldt = as.integer(n), wptr = as.integer(rep(0, 
            d * m)), info = as.integer(0), ptab = as.integer(rep(0, 
            nterms * d)), ldptab = as.integer(nterms), PACKAGE = "fields")
    temp2 <- matrix(temp$tmatrix, nrow = n)
    attr(temp2, "ptab") <- matrix(temp$ptab, nrow = nterms, ncol = d)
    temp2
}
