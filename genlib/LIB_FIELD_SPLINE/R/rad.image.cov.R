"Rad.image.cov" <-
function (ind1, ind2, Y, cov.obj = NULL, setup = FALSE, grid, ...) 
{
    if (is.null(cov.obj)) {
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]
        m <- length(grid$x)
        n <- length(grid$y)
        M <- ceiling2(2 * m)
        N <- ceiling2(2 * n)
        xg <- make.surface.grid(list((1:M) * dx, (1:N) * dy))
        center <- matrix(c((dx * M)/2, (dy * N)/2), nrow = 1, 
            ncol = 2)
        out <- Rad.cov(xg, center, ...)
        out <- as.surface(xg, c(out))$z
        temp <- matrix(0, nrow = M, ncol = N)
        temp[M/2, N/2] <- 1
        wght <- fft(out)/(fft(temp) * M * N)
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
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[ind1])
    }
}
