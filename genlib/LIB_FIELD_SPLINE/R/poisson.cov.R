"poisson.cov" <-
function (x1, x2, eta = 0.2) 
{
    if (missing(x2)) 
        x2 <- x1
#
# dot products of direction cosines  
#
    PP1 <- cbind(cos((x1[, 2] * pi)/180) * cos((x1[, 1] * 
        pi)/180), cos((x1[, 2] * pi)/180) * sin((x1[, 1] * 
        pi)/180), sin((x1[, 2] * pi)/180))
    PP2 <- cbind(cos((x2[, 2] * pi)/180) * cos((x2[, 1] * 
        pi)/180), cos((x2[, 2] * pi)/180) * sin((x2[, 1] * 
        pi)/180), sin((x2[, 2] * pi)/180))
   D  <- (PP1 %*% t(PP2))
con <- (1 - eta^2)^(1.5)
    con/(1 - 2 * eta * D + eta^2)^(1.5)
}
