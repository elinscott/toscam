"Krig.fgcv" <-
function(lam, obj)
{
	#
	# GCV that is leave-one-group out
	#
	lD <- obj$matrices$D * lam
	RSS <- sum(((obj$matrices$u * lD)/(1 + lD))^2)
	MSE <- RSS/length(lD)
	if((obj$N - length(lD)) > 0) {
		MSE <- MSE + obj$pure.ss/(obj$N - length(lD))
	}
	trA <- sum(1/(1 + lD))
	den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(
		lD))
	# If the denominator is negative then flag this as a bogus case
	# by making the GCV function "infinity"
	#
	ifelse(den > 0, MSE/den^2, NA)
}
