"Krig.fdf" <-
function(llam, info)
{
	sum(1/(1 + exp(llam) * info$D)) - info$df
}
