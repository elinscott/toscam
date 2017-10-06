"sreg.trace" <-
function(h, info)
{
	N <- length(info$x)
	#	h <- log(lam)
	temp <- (.Fortran("css",
		h = as.double(h),
		npoint = as.integer(N),
		x = as.double(info$x),
		y = as.double(rep(0, N)),
		wt = as.double(1/sqrt(info$wt)),
		sy = as.double(rep(0, N)),
		trace = as.double(0),
		diag = as.double(rep(0, N)),
		cv = as.double(0),
		ngrid = as.integer(0),
		xg = as.double(0),
		yg = as.double(0),
		job = as.integer(c(3, 0, 0)),
		ideriv = as.integer(0),
		ierr = as.integer(0),PACKAGE="fields")$trace)
	return(temp)
}
