"golden.section.search" <-
function (ax, bx, cx, f, niter = 25, f.extra = NA, tol = 1e-05) 
{
    r <- 0.61803399
    con <- 1 - r
    x0 <- ax
    x3 <- cx
    if (abs(cx - bx) > abs(bx - ax)) {
        x1 <- bx
        x2 <- bx + con * (bx - ax)
    }
    else {
        x2 <- bx
        x1 <- bx - con * (bx - ax)
    }
    f1 <- f(x1, f.extra)

    f2 <- f(x2, f.extra)
    iter <- niter
    for (k in 1:niter) {
#cat( x1,f1, x2,f2, fill=TRUE)
        if (f2 < f1) {
            x0 <- x1
            x1 <- x2
            x2 <- r * x1 + con * x3
            f0 <- f1
            f1 <- f2
            f2 <- f(x2, f.extra)
        }
        else {
            x3 <- x2
            x2 <- x1
            x1 <- r * x2 + con * x0
            f3 <- f2
            f2 <- f1
            f1 <- f(x1, f.extra)
        }
        if (abs(f2 - f1) < tol) {
            iter <- k
            break
        }
    }
    if (f1 < f2) {
        golden <- f1
        xmin <- x1
    }
    else {
        golden <- f2
        xmin <- x2
    }
    list(x = xmin, fmin = golden, iter = iter)
}
