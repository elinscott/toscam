"Krig.fgcv.model" <-
function(lam, obj)
{
	lD <- obj$matrices$D * lam
	MSE <- sum(((obj$matrices$u * lD)/(1 + lD))^2)/length(lD)
	trA <- sum(1/(1 + lD))
	den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(
		lD))
	ifelse(den > 0, obj$shat.pure.error^2 + MSE/den^2, NA)
}
