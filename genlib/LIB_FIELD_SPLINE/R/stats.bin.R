"stats.bin" <-
function (x, y, N = 10, breaks = NULL) 
{
    out <- list()
    if (is.null(breaks)) {
        breaks <- pretty(x, N)
    }
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    test <- describe()
    obj <- matrix(NA, ncol = NBIN, nrow = length(test))
    dimnames(obj) <- list(test, format(1:NBIN))
    obj[, 1] <- describe(y[x <= breaks[2] & x >= breaks[1]])
    for (k in 2:NBIN) {
        obj[, k] <- describe(y[x <= breaks[k + 1] & x > breaks[k]])
    }
    list(centers = centers, breaks = breaks, stats = obj)
}
