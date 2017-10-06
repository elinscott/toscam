"transformx" <-
function (x, scale.type = "unit.sd", x.center, x.scale) 
{
    if (scale.type == "unscaled") {
        x.center <- rep(0, ncol(x))
        x.scale <- rep(1, ncol(x))
    }
    else if (scale.type == "unit.sd") {
        x.center <- apply(x, 2, mean)
        x.scale <- sqrt(apply(x, 2, var))
        x <- scale(x)
    }
    else if (scale.type == "range") {
        x.center <- apply(x, 2, min)
        x.scale <- apply(x, 2, max) - apply(x, 2, min)
        x <- scale(x, center = x.center, scale = x.scale)
    }
    else if (scale.type == "user") {
        if (missing(x.center)) 
            x.center <- apply(x, 2, mean)
        if (missing(x.scale) || length(x.scale) != ncol(x)) 
            stop("Error: x.scale must be a vector of length d")
        x <- scale(x, center = x.center, scale = x.scale)
    }
    else stop(paste("Error: scale.type must be one of", "unit.sd, range, user, unscaled"))
    attr(x, "x.center") <- x.center
    attr(x, "x.scale") <- x.scale
    attr(x, "x.scale.type") <- scale.type
    x
}
