"sreg.df.to.lambda" <-
function(df, x, wt, guess = 1, tol = 1.0000000000000001e-05)
{
	if(is.na(df))
		return(NA)
	n <- length(unique(x))
	info <- list(x = x, wt = wt, df = df)
	if(df > n) {
		warning(" df too large to match a lambda value")
		return(NA)
	}
	h1 <- log(guess)
	########## find upper lambda
	for(k in 1:25) {
		tr <- sreg.trace(h1, info)
		if(tr <= df)
			break
		h1 <- h1 + 1.5
	}
	########## find lower lambda
	##########
	h2 <- log(guess)
	for(k in 1:25) {
		tr <- sreg.trace(h2, info)
		if(tr >= df)
			break
		h2 <- h2 - 1.5
	}
	out <- bisection.search(h1, h2, sreg.fdf, tol = tol, f.extra = info)$
		x
	exp(out)
}
