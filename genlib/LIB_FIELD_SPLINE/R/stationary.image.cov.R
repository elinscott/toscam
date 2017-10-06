stationary.image.cov<- function(
     ind1, ind2, Y, cov.obj = NULL, setup = FALSE, grid, M=NULL, N=NULL,
     Covariance="Matern", Distance="rdist", ...)
{
#
# if cov object is missing then create
# basically need to enlarge domain and find the FFT of the 
# covariance
#
    if (is.null(cov.obj)) {
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]
        m <- length(grid$x)
        n <- length(grid$y)
#
# determine size of padding 
# default is twice domain and will then yeild exact results
# 
        if(is.null(M)) M <- ceiling2(2 * m)
        if(is.null(N)) N <- ceiling2(2 * n)
        xg <- make.surface.grid(list((1:M) * dx, (1:N) * dy))
        center <- matrix(c((dx * M)/2, (dy * N)/2), nrow = 1,
            ncol = 2)
#
# here is where the actual covarinace form is used
# note passed arguments from call for parameters etc.
#
        out <- stationary.cov(xg, center, Covariance=Covariance,
                       Distance=Distance,...)
# coerce to a matrix (image)  
        out <- matrix( c(out),nrow=M, ncol=N)
        temp <- matrix(0, nrow = M, ncol = N)
#
# a simple way to normalize. This could be avoided by 
# translating image from the center ...
#
        temp[M/2, N/2] <- 1
        wght <- fft(out)/(fft(temp) * M * N)
#
# wght is the discrete FFT for the covariance suitable for fast
# multiplication by convolution. 
#
        cov.obj <- list(m = m, n = n, grid = grid, N = N, M = M,
            wght = wght, call = match.call())

        if (setup) {
            return(cov.obj)
        }
    }

    temp <- matrix(0, nrow = cov.obj$M, ncol = cov.obj$N)
    if (missing(ind1)) {
        temp[1:cov.obj$m, 1:cov.obj$n] <- Y
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[1:cov.obj$m,
            1:cov.obj$n])
    }
    else {
        if (missing(ind2)) {
            temp[ind1] <- Y
        }
        else {
            temp[ind2] <- Y
        }

#
# as promised this is a single clean step
#
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[ind1])

    }
}
