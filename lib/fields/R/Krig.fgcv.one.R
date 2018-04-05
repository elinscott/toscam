"Krig.fgcv.one" <-
function(lam, obj)
{
	lD <- obj$matrices$D * lam
	RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
	trA <- sum(1/(1 + lD))
	den <- 1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/obj$N
	# If the denominator is negative then flag this as a bogus case
	# by making the GCV function "infinity"
	#
	ifelse(den > 0, (RSS/obj$N)/den^2, 1e+20)
}
