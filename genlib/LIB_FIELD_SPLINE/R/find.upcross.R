"find.upcross" <-
function (fun, fun.info, upcross.level = 0, guess = 1, tol = 1e-05) 
{
    l1 <- guess
    tr <- 0
    for (k in 1:50) {
        tr <- fun(l1, fun.info) - upcross.level
        if (tr >= 0) 
            break
        else {
            guess <- l1
        }
        l1 <- l1 * 2
    }
    if (tr < 0) {
        warning("Failed to find the upcrossing")
        return(NA)
    }
    tr <- 0
    l2 <- guess
    for (k in 1:50) {
        tr <- fun(l2, fun.info) - upcross.level
        if (tr <= 0) 
            break
        l2 <- l2/2
    }
    if (tr > 0) {
        warning("Failed to find the upcrossing")
        return(NA)
    }
    out <- bisection.search(l2, l1, fun, tol = tol, f.extra = fun.info, 
        upcross.level = upcross.level)$x
    (out)
}
