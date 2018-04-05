"summary.gcv.Krig" <-
function (object, lambda, cost = 1, verbose = FALSE, offset = 0, 
    y = NULL, ...) 
{
    out <- object
    nt <- out$nt
    np <- out$np
    N <- out$N
    D <- out$matrices$D
    if (is.null(y)) {
        u <- out$matrices$u
        shat.pure.error <- out$shat.pure.error
        pure.ss <- out$pure.ss
    }
    else {
        out2 <- Krig.coef(out, y)
        u <- out2$u
        shat.pure.error <- out2$shat.pure.error
        pure.ss <- out2$pure.ss
    }
    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, shat.pure.error = shat.pure.error, 
        offset = offset)
    if (verbose) {
        print(info)
    }
    lambda.est <- rep(NA, 6)
    names(lambda.est) <- c("lambda", "trA", "GCV", "GCV.one", 
        "GCV.model", "shat")
    lambda.est[1] <- lambda
    lambda.est[2] <- Krig.ftrace(lambda, D)
    lambda.est[3] <- Krig.fgcv(lambda, info)
    lambda.est[4] <- Krig.fgcv.one(lambda, info)

    if (!is.na(shat.pure.error)) {
        lambda.est[5] <- Krig.fgcv.model(lambda, info)}

    lambda.est[6] <- sqrt(Krig.fs2hat(lambda, info))
    lambda.est
}

