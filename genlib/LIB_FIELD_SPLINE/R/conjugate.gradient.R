"conjugate.gradient" <-
function (b, multAx, start = NULL, tol = 1e-05, kmax = 25, verbose = TRUE, 
    ...) 
{
    call <- match.call()
    if (is.null(start)) {
        x <- rep(0, length(y))
        r <- b
    }
    else {
        x <- start
        r <- b - multAx(x, ...)
    }
    p <- r
    rho <- rep(NA, kmax + 1)
    rho[1 + (0)] <- sum(r^2)
    test <- sqrt(sum(b^2)) * tol
    niter <- 1
    k <- 0
    for (k in 1:kmax) {
        niter <- niter + 1
        if (k != 1) {
            beta <- rho[1 + (k - 1)]/rho[1 + (k - 2)]
            p <- r + beta * p
        }
        w <- multAx(p, ...)
        alpha <- rho[1 + (k - 1)]/(sum(p * w))
        x <- x + alpha * p
        r <- r - alpha * w
        rho[1 + k] <- sum(r^2)
        if (verbose) {
            cat("iter", k, " crit : ", signif(sqrt(rho[1 + k]), 
                4), signif(test, 4), fill = TRUE)
        }
        if (sqrt(rho[1 + k]) < test) {
            niter <- k + 1
            break
        }
    }
    list(call = call, x = x, residuals = r, niter, conv = list(rho = rho, 
        test = test, maxiter = kmax, niter = niter))
}
