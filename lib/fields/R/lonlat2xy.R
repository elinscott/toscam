"lonlat2xy" <-
function (lnlt, miles = FALSE) 
{
    n <- nrow(lnlt)
    ctr <- c(mean(lnlt[, 1]), mean(lnlt[, 2]))
    ctx <- cbind(lnlt[, 1], rep(ctr[2], n))
    cty <- cbind(rep(ctr[1], n), lnlt[, 2])
    dstx <- c(rdist.earth(matrix(ctr, 1, 2), ctx, miles = FALSE))
    dsty <- c(rdist.earth(matrix(ctr, 1, 2), cty, miles = FALSE))
    cbind(dstx * (2 * (lnlt[, 1] > ctr[1]) - 1), dsty * (2 * 
        (lnlt[, 2] > ctr[2]) - 1))
}
