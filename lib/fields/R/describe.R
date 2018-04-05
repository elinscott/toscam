"describe" <-
function (x) 
{
    lab <- c("N", "mean", "Std.Dev.", "min", "Q1", "median", 
        "Q3", "max", "missing values")
    if (missing(x)) {
        return(lab)
    }
    temp <- rep(0, length(lab))
    xt <- x[!is.na(x)]
    ix <- order(xt)
    n <- length(xt)
    if (!is.numeric(xt) || all(is.na(x))) {
        return(c(n, rep(NA, length(lab) - 2), length(x) - length(xt)))
    }
    if (n == 1) {
        return(c(n, xt[1], NA, rep(xt[1], 5), length(x) - length(xt)))
    }
    else {
        return(c(n, mean(xt), sqrt(var(xt)), min(xt), quantile(xt, 
            c(0.25, 0.5, 0.75)), max(xt), length(x) - length(xt)))
    }
}
