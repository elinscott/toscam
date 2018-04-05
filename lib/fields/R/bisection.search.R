"bisection.search" <-
function (x1, x2, f, tol = 1e-07, niter = 25, f.extra = NA, upcross.level = 0) 
{
    f1 <- f(x1, f.extra) - upcross.level
    f2 <- f(x2, f.extra) - upcross.level
    if (f1 > f2) 
        stop(" f1 must be < f2 ")
    iter <- niter
    for (k in 1:niter) {
        xm <- (x1 + x2)/2
        fm <- f(xm, f.extra) - upcross.level
        if (fm < 0) {
            x1 <- xm
            f1 <- fm
        }
        else {
            x2 <- xm
            f2 <- fm
        }
        if (abs(fm) < tol) {
            iter <- k
            break
        }
    }
    xm <- (x1 + x2)/2
    fm <- f(xm, f.extra) - upcross.level
    list(x = xm, fm = fm, iter = iter)
}
