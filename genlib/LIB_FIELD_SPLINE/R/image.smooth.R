"image.smooth" <-
function (x, wght = NULL, dx = 1, dy = 1, Nwidth = nrow(Y), Mwidth = 
ncol(Y), 
    kernel.function = double.exp, theta = 1, grid = NULL, tol = 1e-08,
     xwidth = NULL, ywidth = NULL, 
    weights = NULL,...) 
{
Y<- x # hack S3
    if (!is.matrix(Y)) {
        stop("Requires a matrix")
    }
    m <- nrow(Y)
    n <- ncol(Y)
    if (!is.null(grid)) {
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]
    }
    if (!is.null(xwidth)) {
        Mwidth <- round(xwidth/dx)
    }
    if (!is.null(ywidth)) {
        Nwidth <- round(ywidth/dy)
    }
    if (is.null(wght)) {
        wght <- image.smooth.setup(nrow = m, ncol = n, Mwidth = Mwidth, 
            Nwidth = Nwidth, dx = dx, dy = dy, kernel.function = kernel.function, 
            theta = theta)
    }
    M <- nrow(wght)
    N <- ncol(wght)
    temp <- matrix(0, nrow = M, ncol = N)
    temp2 <- matrix(0, nrow = M, ncol = N)
    if (!is.null(weights)) {
        temp[1:m, 1:n] <- Y * weights
        temp[is.na(temp)] <- 0
        temp2[1:m, 1:n] <- ifelse(!is.na(Y), weights, 0)
    }
    else {
        temp[1:m, 1:n] <- Y
        temp[is.na(temp)] <- 0
        temp2[1:m, 1:n] <- ifelse(!is.na(Y), 1, 0)
    }
    temp <- Re(fft(fft(temp) * wght, inverse = TRUE))[1:m, 1:n]
    temp2 <- Re(fft(fft(temp2) * wght, inverse = TRUE))[1:m,1:n]
    temp <- ifelse((temp2 > tol), (temp/temp2), NA)
    if (!is.null(grid)) {
        list(x = grid$x, y = grid$y, z = temp)
    }
    else {
        temp
    }
}
