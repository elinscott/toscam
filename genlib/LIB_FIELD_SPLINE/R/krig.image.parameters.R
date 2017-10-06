"krig.image.parameters" <-
function (out)
{
    nt <- out$qr.T$rank
    n <- length(out$yM)
    delta <- qr.qy(out$qr.T, c(rep(0, nt), out$omega2))
    temp <- out$yM - (out$cov.function(out$indexM, , delta, cov.obj = out$cov.obj) +
        out$lambda * out$weightsM * delta)
    beta <- qr.coef(out$qr.T, temp)
    rhohat <- sum(delta * out$yM)/(n - nt)
    rho <- rhohat
    sigma2 <- rho * out$lambda
    shat.MLE <- sqrt(rhohat * out$lambda)
    return(list( beta=beta, delta=delta, rhohat=rhohat, rho=rho, sigma2=sigma2, shat.MLE=shat.MLE))
}

