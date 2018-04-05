"gcv.Krig" <-
function (out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, 
    verbose = FALSE, tol = 1e-05, offset = 0, y = NULL, give.warnings = TRUE, 
    give.warnings.REML = FALSE) 
{
    nt <- out$nt
    np <- out$np
    N <- out$N
    D <- out$matrices$D

# Yet another monster function called by Krig
#    but there just many simple steps ...
#

# if a y data vector is not supplied then 
# use the one in the Krig object 
    if (is.null(y)) {
        u <- out$matrices$u
        shat.pure.error <- out$shat.pure.error
        pure.ss <- out$pure.ss
    }
    else {

#with new data need to update some statistics. 
        out2 <- Krig.coef(out, y=y)
        u <- out2$u
        shat.pure.error <- out2$shat.pure.error
        pure.ss <- out2$pure.ss
    }

    if( verbose) { 
           cat("u used:", fill=TRUE)
           print( u)}

#
# generate a reasonable grid of lambda based on equally spaced
# effective degrees of freedom 

    if (is.na(lambda.grid[1])) {
        temp.df <- seq(nt, (np - offset) * 0.95, , nstep.cv)
        temp.df[1] <- temp.df[1] + 0.001
        for (k in 1:nstep.cv) {
            lambda.grid[k] <- Krig.df.to.lambda(temp.df[k], D)
        }
    }

# make sure that the grid is ordered. 
    lambda.grid <- sort(lambda.grid)

    nl <- length(lambda.grid)
    nd <- length(D)
    V <- V.model <- V.one <- lplike <- trA <- shat <- rep(NA, 
        nl)
    Dl <- rep(NA, nd)
#
# this is small little list used to pass information to the 
# objective functions

    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, shat.pure.error = shat.pure.error, 
        offset = offset)
#
# loop over lambda values for the grid search
 
   for (k in 1:nl) {
#
#       all the wonderful thigns calculated for each lambda
# note the use of the info list. 
        V[k] <- Krig.fgcv(lambda.grid[k], info)
        V.one[k] <- Krig.fgcv.one(lambda.grid[k], info)
        V.model[k] <- Krig.fgcv.model(lambda.grid[k], info)
        lplike[k] <- Krig.flplike(lambda.grid[k], info)
        shat[k] <- sqrt(Krig.fs2hat(lambda.grid[k], info))
        trA[k] <- Krig.ftrace(lambda.grid[k], D)
    }
#
# clean  up as amatrix with all these values. 
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat, 
        lplike)
    gcv.grid <- as.data.frame(gcv.grid)
    names(gcv.grid) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model", 
        "shat", "-Log Profile ")

# at this point the lazy data analyst could stop
# and just pick the minimum over the grid
# but we will now optimze the search producing refined optima
#
 
   if (verbose) 
        print(info)
    if (verbose) 
        print(gcv.grid)

#
# setup output matrix for refined values

    lambda.est <- matrix(ncol = 4, nrow = 6, dimnames = list(c("GCV", 
        "GCV.model", "GCV.one", "RMSE", "pure error", "-Log Profile (REML)"), 
        c("lambda", "trA", "GCV", "shat")))

#
# now step through the many different ways to find lambda

# This is standard GCV w/o replicates, but with replicates this 
# is leave-one-group out. 
#
    lambda.est[1, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV, 
        Krig.fgcv, tol = tol, verbose = verbose, give.warnings = give.warnings)
#
# If replicates then try GCV related to replicate group means
# 
   if (!is.na(shat.pure.error)) {
        temp <- gcv.grid$GCV.model
        lambda.est[2, 1] <- Krig.find.gcvmin(info, lambda.grid, 
            temp, Krig.fgcv.model, tol = tol, verbose = verbose, 
            give.warnings = give.warnings)
    }
#
# Another version of GCV this is really leave an observation out

    lambda.est[3, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV.one, 
        Krig.fgcv.one, tol = tol, verbose = verbose, give.warnings = give.warnings)
#
# REML 
    lambda.est[6, 1] <- Krig.find.REML(info, lambda.grid, gcv.grid$"-Log Profile", 
        Krig.flplike, tol = tol, verbose = verbose, give.warnings = give.warnings.REML)
    if (verbose) {
        cat(" mle estimate", lambda.est[6, 1], fill = TRUE)
    }
# 
# This is the most obscure Perry Haaland likes this one ....
# Find the values of lambda so that the 
# usual estimate of sigma**2 is equal to either the 
# value estimated from replicates and  from an external value (rmse) if
# this has been passed. 
# in some sense this a moment matching choice for lambda
    lambda.rmse <- NA
    lambda.pure.error <- NA

# matching external RMSE
    if (!is.na(rmse)) {
        if (all(gcv.grid$shat < rmse) | all(gcv.grid$shat > rmse)) {
            guess <- NA
        }
        else {
            guess <- max(gcv.grid$lambda[gcv.grid$shat < rmse])
        }

        if (verbose) {
            print(rmse)
            print(guess)
        }
        if (!is.na(guess)) {
            lambda.rmse <- find.upcross(Krig.fs2hat, info, upcross.level = rmse^2, 
                guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
        }
        else {
            if (give.warnings) {
                warning("Value of rmse is outside possible range")
            }
        }
    }
#
# matching estimate of sigma from reps. 

    if (!is.na(shat.pure.error)) {
        if (all(gcv.grid$shat < shat.pure.error) | all(gcv.grid$shat > 
            shat.pure.error)) {
            guess <- NA
        }
        else {
            guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])
        }

        if (!is.na(guess) & (guess != -Inf)) {
            lambda.pure.error <- find.upcross(Krig.fs2hat, info, 
                upcross.level = shat.pure.error^2, guess = guess, 
                tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
        }
        else {
            if (give.warnings) {
                warning("Value of pure error estimate  is outside possible range")
            }
        }
    }

    if (verbose) {
        cat(" Done with refining", fill = TRUE)
    }

#
# OK done with all six methods. 
# fill in return matrix with all the right stuff 


    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            lambda.est[k, 2] <- Krig.ftrace(lam, D)
            if (k == 1 | k > 3) {
                lambda.est[k, 3] <- Krig.fgcv(lam, info)
            }
            if (k == 2) {
                lambda.est[k, 3] <- Krig.fgcv.model(lam, info)
            }
            if (k == 3) {
                lambda.est[k, 3] <- Krig.fgcv.one(lam, info)
            }
            lambda.est[k, 4] <- sqrt(Krig.fs2hat(lam, info))
        }
              if( verbose){
                  cat( "lam.est ",k, fill=TRUE)
                  cat( lambda.est[k,], fill=TRUE)}
    }
    lam.ml <- lambda.est[6, 1]
    lambda.est[6, 2] <- Krig.ftrace(lam.ml, D)
    lambda.est[6, 3] <- Krig.fgcv(lam.ml, info)
    lambda.est[6, 4] <- sqrt(Krig.fs2hat(lam.ml, info))

# Note that the preferred estimate by default is 
# leave-one (obs or group)-out GCV 

  list(gcv.grid = gcv.grid, 
      lambda.est = lambda.est, lambda.best = lambda.est[1, 1])
}

