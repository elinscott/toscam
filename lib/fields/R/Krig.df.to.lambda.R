"Krig.df.to.lambda" <-
function (df, D, guess = 1, tol = 1e-05) 
{
    if (is.list(D)) {
        D <- D$matrices$D
    }
    if (is.na(df)) 
        return(NA)
    if (df < sum(D == 0)) {
        warning("df too small to match with a lambda value")
        return(NA)
    }
    if (df > length(D)) {
        warning(" df too large to match a lambda value")
        return(NA)
    }
    l1 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l1 * D))
        if (tr <= df) 
            break
        l1 <- l1 * 4
    }
    l2 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l2 * D))
        if (tr >= df) 
            break
        l2 <- l2/4
    }
    info <- list(D = D, df = df, N = length(D))
    out <- bisection.search(log(l1), log(l2), Krig.fdf, tol = tol, 
        f.extra = info)$x
    +exp(out)
}

