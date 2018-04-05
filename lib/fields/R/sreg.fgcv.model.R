"sreg.fgcv.model" <-
function(lam, obj)
{
	sreg.fit(lam, obj)$gcv.model
}
