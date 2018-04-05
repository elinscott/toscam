"print.summary.krig.image" <-
function (x, ...) 
{
    digits <- x$digits
    c1 <- "Number of Observations:"
    c2 <- x$num.observation
    c1 <- c(c1, "Degree of polynomial null space ( base model):")
    c2 <- c(c2, x$m - 1)
    c1 <- c(c1, "MLE sigma ")
    c2 <- c(c2, format(signif(x$shat.MLE, digits)))
    c1 <- c(c1, "MLE rho ")
    c2 <- c(c2, format(signif(x$rhohat, digits)))
    c1 <- c(c1, "Scale used for covariance (rho)")
    c2 <- c(c2, signif(x$rho, digits))
    c1 <- c(c1, "Scale used for nugget (sigma^2)")
    c2 <- c(c2, signif(x$sigma2, digits))
    c1 <- c(c1, "lambda (sigma2/rho)")
    c2 <- c(c2, signif(x$lambda, digits))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    res.quantile <- x$res.quantile
    names(res.quantile) <- c("min", "1st Q", "median", "3rd Q", 
        "max")
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    if ((x$correlation.model)) {
        cat("Y is standardized before spatial estimate is found", 
            fill = TRUE)
    }
    cat(" Residuals:", "\n")
    print(signif(res.quantile, digits))
    cat("Covariance function name:", x$cov.function, fill = TRUE)
    invisible(x)
}
