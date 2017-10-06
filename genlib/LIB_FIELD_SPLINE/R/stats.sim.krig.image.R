"stats.sim.krig.image" <-
function (obj) 
{
    n <- length(obj$grid$y)
    m <- length(obj$grid$x)
    out1 <- out2 <- matrix(0, m, n)
    out3 <- out1[-1, ]
    out4 <- t(out3)
    k <- 1
    temp <- out1
    N <- length(obj$out)
    for (k in 1:N) {
        temp <- obj$out[[k]]
        out1 <- out1 + temp
        out2 <- temp^2 + out2
        out3 <- temp[1:(m - 1), ] * temp[2:m, ] + out3
        out4 <- temp[, 1:(n - 1)] * temp[, 2:n] + out4
        NULL
    }
    Hcor <- out3/sqrt(out2[1:(m - 1), ] * out2[2:m, ])
    Vcor <- out4/sqrt(out2[, 1:(n - 1)] * out2[, 2:n])
    out1 <- out1/N
    return(grid = obj$grid, var = (out2 - N * out1^2)/(N - 1), 
        Vcor, Hcor, Vcor, mean = out1)
}
