 "smooth.2d" <-
function (Y, ind = NULL, weight.obj = NULL, setup = FALSE, grid = NULL, 
    x = NULL, nrow = 64, ncol = 64, surface = TRUE, cov.function = 
gauss.cov, 
    Mwidth = NULL, Nwidth = NULL, ...) 
{
    temp <- as.image(Y, ind, grid = grid, nrow = nrow, ncol = ncol, 
        x = x)
    Y <- temp$z
    NN <- temp$N
    grid <- list(x = temp$x, y = temp$y)
    if (is.null(weight.obj)) {
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]
        m <- length(grid$x)
        n <- length(grid$y)
        if (is.null(Mwidth)) 
            M <- 2 * m
        else {
            M <- m + Mwidth
        }
        if (is.null(Nwidth)) 
            N <- 2 * n
        else {
            N <- n + Nwidth
        }
        xg <- make.surface.grid(list((1:M) * dx, (1:N) * dy))
        center <- matrix(c((dx * M)/2, (dy * N)/2), nrow = 1, 
            ncol = 2)
        out <- cov.function(xg, center, ...)
        out <- as.surface(xg, c(out))$z
        temp <- matrix(0, nrow = M, ncol = N)
        temp[M/2, N/2] <- 1
        wght <- fft(out)/(fft(temp) * M * N)
        weight.obj <- list(m = m, n = n, grid = grid, N = N, 
            M = M, wght = wght, call = match.call())
        if (setup) {
            return(weight.obj)
        }
    }
    temp <- matrix(0, nrow = weight.obj$M, ncol = weight.obj$N)
    temp[1:m, 1:n] <- Y
    temp[is.na(temp)] <- 0
    temp2 <- Re(fft(fft(temp) * weight.obj$wght, inverse = TRUE))[1:weight.obj$m, 
        1:weight.obj$n]
    temp <- matrix(0, nrow = weight.obj$M, ncol = weight.obj$N)
    temp[1:m, 1:n] <- NN
    temp[is.na(temp)] <- 0
    temp3 <- Re(fft(fft(temp) * weight.obj$wght, inverse = TRUE))[1:weight.obj$m, 
        1:weight.obj$n]
    if (!surface) 
        (temp2/temp3)
    else {
        list(x = weight.obj$grid$x, y = weight.obj$grid$y, z = temp2/temp3, 
            index = ind)
    }
}
