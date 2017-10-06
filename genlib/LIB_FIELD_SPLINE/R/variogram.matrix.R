"variogram.matrix" <-
function (dat, R = 5,dx) 
{
    shift.fct <- function(d) {
        shift.x <- c(0:d)
        shift.y <- sqrt(round(d^2) - shift.x^2)
        indices <- ((as.integer(shift.y))^2 == (round(d^2) - 
            shift.x^2))
        shift.y <- shift.y[indices]
        shift.x <- shift.x[indices]
        return(list(x=shift.x, y=shift.y))
    }
    n <- ncol(dat)
    m <- nrow(dat)
    a <- as.integer(R/sqrt(2)) + 1
    if (a > max(m, n)) {
        cat("\n a=", a, " has be larger than max(m,n) \n")
        return()
    }
    ind <- cbind(rep(1:a, a), rep(1:a, rep(a, a)))
    x.diff <- outer(ind[, 1], ind[, 1], "-")
    y.diff <- outer(ind[, 2], ind[, 2], "-")
    d.matrix <- sqrt(x.diff^2 + y.diff^2)
    d <- d.matrix[, 1]
    indices <- c(0)
    i <- 0
    l.d <- length(d)
    while ((i + 1) * a < l.d) {
        indices <- c(indices, c((i * a + i + 1):((i + 1) * a)))
        i <- i + 1
    }
    indices <- c(indices, (i * a + i + 1):l.d)
    indices <- indices[-2]
    l.i <- length(indices)
    d <- d[indices]
    d <- sort(d)
    vgram <- rep(NA, length(d))
    for (i in 1:length(d)) {
        sum.temp <- 0
        l.temp <- 0
        shift <- shift.fct(d[i])
        for (j in 1:length(shift$x)) {
            h.x <- shift$x[j]
            h.y <- shift$y[j]
            sum.temp <- sum.temp + sum((dat[1:(m - h.x), 1:(n - 
                h.y)] - dat[(1:(m - h.x) + h.x), (1:(n - h.y) + 
                h.y)])^2)
            l.temp <- l.temp + ((n - h.x) * (m - h.y))
        }
        vgram[i] <- (0.5 * sum.temp)/l.temp
    }
    list(vgram=vgram, d=d*dx)
}
