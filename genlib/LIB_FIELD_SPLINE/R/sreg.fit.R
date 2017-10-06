"sreg.fit" <-
function(lam, info)
{
	np <- info$np
	N <- info$N
	nt <- 2
	if(is.null(info$cost)) {
		cost <- 1
	}
	else {
		cost <- info$cost
	}
	if(is.null(info$offset)) {
		offset <- 0
	}
	else {
		offset <- info$offset
	}
	if(is.null(info$shat.pure.error)) {
		shat.pure.error <- 0
	}
	else {
		shat.pure.error <- info$shat.pure.error
	}
	if(is.null(info$pure.ss)) {
		pure.ss <- 0
	}
	else {
		pure.ss <- info$pure.ss
	}
	#print(np)
	#	NOTE h <- log(lam)
	temp <- .Fortran("css",
		h = as.double(log(lam)),
		npoint = as.integer(np),
		x = as.double(info$x),
		y = as.double(info$y),
		wt = as.double(sqrt(1/info$wt)),
		sy = as.double(rep(0, np)),
		trace = as.double(0),
		diag = as.double(rep(0, np)),
		cv = as.double(0),
		ngrid = as.integer(0),
		xg = as.double(0),
		yg = as.double(0),
		job = as.integer(c(3, 0, 0)),
		ideriv = as.integer(0),
		ierr = as.integer(0),
                PACKAGE="fields")
	rss <- sum((temp$sy - info$y)^2 * info$wt)
	MSE <- rss/np
	if((N - np) > 0) {
		MSE <- MSE + pure.ss/(N - np)
	}
	trA <- temp$trace
	den <- (1 - (cost * (trA - nt - offset) + nt)/np)
	den1 <- (1 - (cost * (trA - nt - offset) + nt)/N)
	# If the denominator is negative then flag this as a bogus case
	# by making the GCV function "infinity"
	#
	shat <- sqrt((rss + pure.ss)/(N - trA))
	GCV <- ifelse(den > 0, MSE/den^2, NA)
	gcv.model <- ifelse(
(den > 0)&( (N-np)>0), pure.ss/(N - np) + (rss/np)/(den^2), NA)
	gcv.one <- ifelse(den > 0, ((pure.ss + rss)/N)/(den1^2), NA)
	list(trace = trA, gcv = GCV, rss = rss, shat = shat, gcv.model = 
		gcv.model, gcv.one = gcv.one)
}
