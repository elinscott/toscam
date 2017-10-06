"periodic.cov.1d" <-
function (x1, x2, a, b) 
{
    pairdist <- rdist(x1, x2)%%(b - a)/(b - a)
    bernoulli4 <- pairdist^4 - 2 * (pairdist^3) + pairdist^2 - 
        1/30
    cov <- ((-1) * bernoulli4)/(4 * 3 * 2)
    cov <- cov + (max(cov) - min(cov))/4
    cov <- cov/max(abs(cov))
    return(cov)
}
