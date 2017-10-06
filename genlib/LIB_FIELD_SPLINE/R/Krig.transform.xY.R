Krig.transform.xY<- function(obj,knots, verbose=FALSE ){

# find all replcates and  collapse to unique locations and mean response 
# and pooled variances and weights. 

out<- Krig.replicates( obj, verbose=verbose)

#
# save information about knots.

    if (is.na(knots[1])) {
        out$knots <- out$xM
        out$mle.calc <- TRUE
        out$knot.model <- FALSE
    }
    else {
        out$mle.calc <- FALSE
        out$knot.model <- TRUE
        out$knots<- knots
    }

#
# scale x, knot locations and  save transformation info
#
    out$xM <- transformx(out$xM, obj$scale.type, obj$x.center, obj$x.scale)
    out$transform <- attributes(out$xM)
    out$knots <-   scale(out$knots, center = out$transform$x.center,
                   scale = out$transform$x.scale)
#
#
#verbose block
#
    if (verbose) {
        cat("transform", fill = TRUE)
        print(out$transform)
    }
    if (verbose) {
        cat("knots in transformed scale", fill = TRUE)
        print(knots)
    }

return( out)
}

