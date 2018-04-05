"krig.image" <-
function (x, Y, cov.function, m = NULL, n = NULL, lambda = 0, 
    start = NULL, tol = 1e-05, kmax = 25, cov.obj = NULL, grid = NULL, 
    weights = rep(1, length(Y)), verbose = FALSE, conv.verbose = FALSE, 
    expand = 1, ...) 
{
    out <- list()
    out$call <- match.call()
    out$cov.function <- cov.function
    na.ind <- is.na(Y)
    out$na.ind <- na.ind
    out$xraw <- x[!na.ind, ]
    out$y <- Y[!na.ind]
    out$N <- length(out$y)
    out$weights <- weights[!na.ind]
    out$lambda <- lambda
    if (is.null(grid)) {
        if (is.null(cov.obj)) {
            x1r <- range(out$xraw[, 1])
            x2r <- range(out$xraw[, 2])
            if (verbose) {
                print(x1r)
                print(x2r)
            }
            out$grid <- list(x = seq(x1r[1], x1r[2], , m), y = seq(x2r[1], 
                x2r[2], , n))
        }
        else out$grid <- cov.obj$grid
    }
    else out$grid <- grid
    if (verbose) {
        print(out$grid)
    }
    if (is.null(cov.obj)) {
        cov.obj <- cov.function(grid = out$grid, setup = TRUE, ...)
    }
    out$cov.obj <- cov.obj
    out <- c(out, discretize.image(out$xraw, grid = out$grid, 
        expand = expand))
    if (verbose) {
        cat("row and columns of grid", fill = TRUE)
        print(c(out$m, out$n))
    }
    out <- c(out, list(x = cbind(out$grid$x[out$index[, 1]], 
        out$grid$y[out$index[, 2]])))
    if (verbose) {
        cat("results of finding indices X,y,weights, indices", 
            fill = TRUE)
        print(cbind(out$x, out$y, out$weights, out$index))
    }
    out <- c(out, Krig.replicates(out, verbose = verbose))
    out$indexM <- out$index[out$uniquerows, ]
    if (verbose) {
        cat("results of collapsing replicates", fill = TRUE)
        print(cbind(out$xM, out$yM, out$weightsM, out$indexM))
    }
    out$qr.T <- qr(cbind(rep(1, nrow(out$indexM)), out$xM))
    if (is.null(start)) {
        start <- rep(0, length(out$yM) - out$qr.T$rank)
    }
    out$multAx <- function(x, loc, cov.function, cov.obj, D, 
        qr.T) {
        temp <- qr.qy(qr.T, c(rep(0, qr.T$rank), x))
        qr.q2ty(qr.T, cov.function(loc, , temp, cov.obj = cov.obj) + 
            D * temp)
    }
    Q2TY <- qr.q2ty(out$qr.T, out$yM)
    cat(" iterative solution of linear system", fill = TRUE)
    out2 <- conjugate.gradient(b = Q2TY, out$multAx, start = start, 
        tol = tol, kmax = kmax, D = out$lambda * out$weightsM, 
        loc = out$indexM, cov.function = out$cov.function, cov.obj = out$cov.obj, 
        qr.T = out$qr.T, verbose = conv.verbose)
    if (verbose) {
        cat("convergence info from conjugate gradient", fill = TRUE)
        print(out2$conv)
    }
    out$omega2 <- out2$x
    out$converge <- out2$conv
    if (out$converge$niter >= kmax) {
        warning("Exceeded max number of\niterations for congate gradient.")
    }
    out <- c(out, krig.image.parameters(out))
    class(out) <- "krig.image"
    if (verbose) {
        cat("estimated parameters delta, then beta", fill = TRUE)
        print(out$delta)
        print(out$beta)
    }
    if (verbose) {
        cat("estimated parameters sigma2 and rho", fill = TRUE)
        print(out$sigma2)
        print(out$rho)
    }
    xg <- make.surface.grid(out$grid)
    ig <- make.surface.grid(list(1:out$m, 1:out$n))
    z1 <- out$cov.function(ig, out$indexM, Y = out$delta, cov.obj = out$cov.obj)
    z2 <- out$beta[1] + xg %*% out$beta[2:3]
    out$surface <- as.surface(xg, z1 + z2)
    if (verbose) {
        print(out$surface$x)
        print(out$surface$y)
    }
    if (verbose) 
        print(out$rep.info)
    out$fitted.values <- predict(out, out$xraw)
    out$residuals <- out$y - out$fitted.values
    out
}
