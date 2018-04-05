"Krig.find.REML" <-
function (info, lambda.grid, llike, llike.fun, tol, verbose = TRUE,
    give.warnings = FALSE) 
{
#
# NOTE give.warnings set to FALSE to avoid numerous messages for
# the standard fields examples.  
#
    ind <- !is.na(llike)
    lambda.grid <- lambda.grid[ind]
    llike <- llike[ind]
    nstep.cv <- length(lambda.grid)
    il <- order(llike)[1]
    lambda.llike <- lambda.grid[il]
    llike.raw <- min(llike)
    if (verbose) {
        cat("Results of coarse search lambda and  restricted Log Likelihood:", 
lambda.llike[il], 
            llike.raw, fill = TRUE)
    }
    if ((il > 1) & (il < nstep.cv)) {
        out <- golden.section.search(lambda.grid[il - 1], lambda.grid[il], 
            lambda.grid[il + 1], llike.fun, f.extra = info, tol = tol * 
                llike.raw)
        return(out$x)
    }
    else {
        if (give.warnings) {
            warning("Search for REML estimate of smoothing paramter gives a 
maximum at the endpoints of the grid search")
        }
        return(lambda.llike)
    }
}

