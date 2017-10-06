"sim.krig.image" <-
function (out, nreps = 10) 
{
    temp1 <- list()
    temp2 <- as.list(1:nreps)
    class(temp1) <- "sim.krig.image"
    temp1$call <- out$call
    temp1$grid <- out$grid
    hold2 <- look <- out$surface$z
    temp <- krig.image.parameters(out)
    rho <- temp$rho
    shat.MLE <- temp$shat.MLE
    for (k in 1:nreps) {
        look <- sqrt(rho) * sim.rf(out$cov.obj)
        Y <- c(look[out$index]) + rnorm(out$N) * shat.MLE
        hold2 <- krig.image(out$x, Y, lambda = out$lambda, cov.function = out$cov.function, 
            cov.obj = out$cov.obj, kmax = 100, start = NULL)$surface$z
        temp2[[k]] <- look - hold2 + out$surface$z
        NULL
    }
    temp1$out <- temp2
    return(temp1)
}
