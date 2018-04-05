"gcv.sreg" <-
function (out, lambda.grid = NA, cost = 1, nstep.cv = 20, rmse = NA, 
    offset = 0, trmin = NA, trmax = NA, verbose = TRUE, tol = 1e-05, 
    find.min = TRUE, method = "GCV") 
{
    shat.pure.error <- out$shat.pure.error
    pure.ss <- out$pure.ss
    nt <- 2
    np <- out$np
    N <- out$N
    out$cost <- cost
    out$offset <- offset

# find good end points for lambda coarse grid. 

    if (is.na(trmin)) 
        trmin <- 2.05
    if (is.na(trmax)) 
        trmax <- out$np * 0.95
    if (is.na(lambda.grid[1])) {
        l2 <- sreg.df.to.lambda(trmax, out$x, out$wt)
        l1 <- sreg.df.to.lambda(trmin, out$x, out$wt)
        lambda.grid <- exp(seq(log(l2), log(l1), , nstep.cv))
    }
    if (verbose) {
        cat( "endpoints of coarse lamdba grid", fill=TRUE)
        cat(l1, l2, fill = TRUE)
    }

# build up table of coarse grid serach results for lambda

    nl <- length(lambda.grid)
    V <- V.model <- V.one <- trA <- MSE <- RSS.model <- rep(NA, 
        nl)

#   loop through lambda's and compute various qaunties related to 
#   lambda and the fitted spline. 
    for (k in 1:nl) {
        temp <- sreg.fit(lambda.grid[k], out)
        RSS.model[k] <- temp$rss
        trA[k] <- temp$trace
        V[k] <- temp$gcv
        V.one[k] <- temp$gcv.one
        V.model[k] <- temp$gcv.model
    }
    RSS <- RSS.model + pure.ss
    shat <- sqrt(RSS/(N - trA))
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat)
    dimnames(gcv.grid) <- list(NULL, c("lambda", "trA", "GCV", 
        "GCV.one", "GCV.model", "shat"))

    if (verbose) {
        cat( "Results of coarse grid search", fill=TRUE)
        print(gcv.grid)
    }
    if (!find.min) {
        return(list(gcv.grid = gcv.grid))
    }
    lambda.est <- matrix(ncol = 4, nrow = 5, dimnames = list(c("GCV", 
        "GCV.model", "GCV.one", "RMSE", "pure error"), c("lambda", 
        "trA", "GCV", "shat")))

# now do various refinements for different flavors of finding 
# a good value for lambda the smoothing parameter

##### traditional leave-one-out
    lambda.est[1, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[, 
        "GCV"], sreg.fgcv, tol = tol, verbose = verbose)
    if (verbose) {
        cat("leave one out GCV lambda", fill=TRUE)
        cat(lambda.est[1, 1], fill = TRUE)
    }

##### using GCV criterion adjusting for replicates
    if (!is.na(shat.pure.error)) {
        lambda.est[2, 1] <- Krig.find.gcvmin(out, lambda.grid, 
            gcv.grid[, "GCV.model"], sreg.fgcv.model, tol = tol, 
            verbose = verbose)
        if (verbose) {
            cat("results of GCV.model search", fill=TRUE) 
            cat(lambda.est[2, 1], fill = TRUE)
        }
    }
    lambda.est[3, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[, 
        "GCV.one"], sreg.fgcv.one, tol = tol, verbose = verbose)

##### matching an external value of RMSE
    lambda.rmse <- NA
    lambda.pure.error <- NA
    if (!is.na(rmse)) {
        guess <- max(gcv.grid[gcv.grid[, "shat"] < rmse, "lambda"])

        if (verbose) {
            cat("trying to matching with RMSE", fill=TRUE) 
            cat("rmse", rmse, "guess", guess, fill=TRUE)
        }
   
     if (!is.na(guess)) {
            lambda.rmse <- find.upcross(sreg.fs2hat, out, upcross.level = rmse^2, 
                guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
   
        }
        else {
            warning("Value of rmse is outside possible range")
        }
    }

##### matching  sigma estimated from the replicates. 
    if (!is.na(shat.pure.error)) {
#       all shats smaller than pure error estimate
        guess <- gcv.grid[, "shat"] < shat.pure.error 

#       set to NA a bad guess!
        if( any( guess) ) { 
             guess<- max(  gcv.grid[guess, "lambda"]) }
        else{
           guess<- NA}

        if (verbose) {
            cat("#### trying to matching with sigma from pure error", 
fill=TRUE)
            cat("shat.pure", shat.pure.error, "guess", guess, fill=TRUE)
        }

        if (!is.na(guess)) {
            lambda.pure.error <- find.upcross(sreg.fs2hat, out, 
                upcross.level = shat.pure.error^2, guess = guess, 
                tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
            cat("results of matching with pure error sigma", fill=TRUE) 
        }
        else {
            warning("Value of pure error estimate   
                             is outside possible range")
        }
    }
    if (verbose) {
        cat("All forms of estimated lambdas so far", fill=TRUE)
        print(lambda.est)
    }
    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            temp <- sreg.fit(lam, out)
            lambda.est[k, 2] <- temp$trace
            if ((k == 1) | (k > 3)) {
                lambda.est[k, 3] <- temp$gcv
            }
            if (k == 2) {
                lambda.est[k, 3] <- temp$gcv.model
            }
            if (k == 3) {
                lambda.est[k, 3] <- temp$gcv.one
            }
            lambda.est[k, 4] <- temp$shat
        }
    }
    list(gcv.grid = gcv.grid, lambda.est = lambda.est, lambda.best = lambda.est[method, 
        1])
}

