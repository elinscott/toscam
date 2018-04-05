"summary.sreg" <-
function(object, digits = 4, ...)
{
        x<- object # hack S3
	if(length(x$lambda) > 1) {
		stop("Can't do a summary on an object with a grid of smoothing\nparameters"
			)
	}
	summary <- list(call = x$call, num.observation = length(x$residuals),
		enp = x$trace, nt = x$nt, res.quantile = quantile(x$residuals,
		seq(0, 1, 0.25)), shat.GCV = x$shat.GCV, m = x$m, lambda = x$
		lambda, cost = x$cost, num.uniq = length(x$y), np = x$np, 
		method = x$method, lambda.est = x$lambda.est[!is.na(x$
		lambda.est[, 1]),  ], shat.pure.error = x$shat.pure.error)
	class(summary) <- "summary.sreg"
	summary$covariance <- cor(x$fitted.values * sqrt(x$wt.raw), (x$yraw) *
		sqrt(x$wt.raw))^2
	hold <- (sum((x$yraw - mean(x$yraw))^2) - sum(x$residuals^2))/(sum(
		(x$yraw - mean(x$yraw))^2))
	summary$adjr2 <- 1 - ((length(x$residuals) - 1)/(length(x$residuals) -
		x$eff.df)) * (1 - hold)
	summary$digits <- digits
	summary$sum.gcv.lambda <- summary.gcv.sreg(x, x$lambda)
	summary
}
