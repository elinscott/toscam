"print.sreg" <-
function(x, ...)
{
	if(length(x$lambda) > 1) {
		c1 <- "Number of Observations:"
		c2 <- (x$N)
		c1 <- c(c1, "Number of values of lambda in grid:")
		c2 <- c(c2, length(x$lambda))
		sum <- cbind(c1, c2)
	}
	else {
		digits <- 4
		N <- x$N
		c1 <- "Number of Observations:"
		c2 <- (x$N)
		c1 <- c(c1, "Effective degrees of freedom:")
		c2 <- c(c2, format(round(x$trace, 1)))
		c1 <- c(c1, "Residual degrees of freedom:")
		c2 <- c(c2, format(round(x$N - x$trace, 1)))
		c1 <- c(c1, "Residual root mean square:")
		c2 <- c(c2, format(signif(sqrt(sum(x$residuals^2)/N), 4)))
		c1 <- c(c1, "Lambda ( smoothing parameter)")
		c2 <- c(c2, format(signif((x$lambda), 4)))
		sum <- cbind(c1, c2)
	}
	dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
	cat("Call:\n")
	dput(x$call)
	print(sum, quote = FALSE)
	invisible(x)
}
