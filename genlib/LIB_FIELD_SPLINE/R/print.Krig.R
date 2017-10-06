"print.Krig" <-
function (x, digits = 4, ...) 
{
    c1 <- "Number of Observations:"
    c2 <- length(x$residuals)
#
# print out null space poly info only if "m" is used
    
    if( !is.null( x$args.null$m) ) {
    c1 <- c(c1, "Degree of polynomial null space ( base model):")
    c2 <- c(c2, x$m - 1) }
   
    c1 <- c(c1, "Number of parameters in the null space")
    c2 <- c(c2, x$nt)
    c1<- c( c1,"Parameters for fixed spatial drift")
    c2<- c( c2, sum(x$ind.drift))
    c1 <- c(c1, "Model degrees of freedom:")
    c2 <- c(c2, format(round(x$eff.df, 1)))
    c1 <- c(c1, "Residual degrees of freedom:")
    c2 <- c(c2, format( round(length(x$residuals) - x$eff.df,1) ) )
    c1 <- c(c1, "GCV estimate for sigma:")
    c2 <- c(c2, format(signif(x$shat.GCV, digits)))
    c1 <- c(c1, "MLE for sigma:")
    c2 <- c(c2, format(signif(x$shat.MLE, digits)))

    c1 <- c(c1, "MLE for rho:")
    c2 <- c(c2, format(signif(x$rhohat, digits)))

    c1 <- c(c1, "lambda")
    c2 <- c(c2, format(signif(x$lambda, 2)))

    c1 <- c(c1, "User rho")
    c2 <- c(c2, format(signif(x$rho, digits)))

    c1 <- c(c1, "User sigma^2")
    c2 <- c(c2, format(signif(x$sigma2, digits)))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    invisible(x)
}

