"vgram.matrix" <-
function (dat, R = 5, nsum = 1:8, collapse = TRUE, dx=1,dy=1) 
{

    if (collapse& (dx==dy)) {
        variogram.matrix(dat, round(R/dx),dx)
    }
    else {
        N <- ncol(dat)
        M <- nrow(dat)
        m <- round(R/dx)
        n <- round(R/dy)

        ind <- cbind(rep(0:m, n + 1), rep(0:n, rep(m + 1, n + 
            1)))
        d <- sqrt((dx*ind[, 1])^2 + (dy*ind[, 2])^2)
        ind <- ind[(d > 0) & (d <= R), ]
        d <- d[(d > 0) & (d <= R)]
        ind <- ind[order(d), ]
        d <- sort(d)
        nbin <- nrow(ind)
        ns <- length(nsum)
        hold <- matrix(NA, ncol = ns, nrow = nbin)
        hold2 <- rep(NA, nbin)
        for (k in 1:nbin) {
            m1 <- M - ind[k, 1]
            m2 <- ind[k, 1] + 1
            n1 <- N - ind[k, 2]
            n2 <- ind[k, 2] + 1
            hold[k, ] <- c(describe(0.5 * (dat[1:m1, 1:n1] - 
                dat[m2:M, n2:N])^2))[nsum]
            hold2[k] <- mean((0.5 * (abs(dat[1:m1, 1:n1] - dat[m2:M, 
                n2:N]))^0.5))
        }
        cst <- (0.457 + 0.494/nbin)
        hold2 <- hold2^4/cst
        list(d = d, ind = ind, stats = hold, vgram = hold[, 2], 
            vgram.robust = hold2)
    }
}
