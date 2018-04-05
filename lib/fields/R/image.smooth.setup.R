"image.smooth.setup" <-
function (x, nrow = 64, ncol = 64, dx = 1, dy = 1, kernel.function = 
double.exp, 
    theta = 1, Mwidth = nrow, Nwidth = ncol, ...) 
{
    m <- nrow
    n <- ncol
    M2 <- round((m + Mwidth)/2)
    N2 <- round((n + Nwidth)/2)
    M <- 2 * M2
    N <- 2 * N2
    xi <- (1:M2) * dx
    xi <- xi/theta
    yi <- (1:N2) * dy
    yi <- yi/theta
    out <- kernel.function((matrix(xi, M2, N2)^2 + matrix(yi, 
        M2, N2, byrow = TRUE)^2), ...)/theta
    out <- cbind(out, out[, N2:1])
    out <- rbind(out, out[M2:1, ])
    fft(out)/(M * N)
}
