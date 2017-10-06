"sreg" <-
function(x, y, lam = NA, df = NA, offset = 0, wt = rep(1, length(x)), cost = 1,
	nstep.cv = 80, find.diagA = TRUE, trmin = 2.01, trmax = 
	length(unique(x)) * 0.95, lammin = NA, lammax = NA,
	verbose = FALSE, do.cv = TRUE, method = "GCV", rmse = NA, lambda = NA)
{
	call <- match.call()
	out <- list()
	out$call <- match.call()
	class(out) <- c("sreg")
	out$xraw <- x
	out$yraw <- y
	out$wt.raw <- wt
	out$cost <- cost
	out$offset <- offset
	out$method <- method
	if(length(x) != length(y))
		stop(" X and Y do not match")
	N <- length(y)
	out$N <- N
	out$nt <- 2
	#
	# figure out if the GCV function should be minimized
	# and that value of lambda used for the estimate
	#
	if(!is.na(lambda[1])) {
		lam <- lambda
	}
	#
	#
	if(is.na(lam[1]) & is.na(df[1])) {
		eval.cv <- TRUE
		find.min <- TRUE
	}
	else {
		find.min <- FALSE
		eval.cv <- FALSE
	}
	#
	## find duplicate rows of the x vector 
	## first find integer tags to indicate replications
	rep.info <- cat.matrix(matrix(x, ncol = 1))
	out$rep.info <- rep.info
	if(verbose) {
		cat("rep.info", fill = TRUE)
		print(rep.info)
	}
	if(max(rep.info) == N) {
		shat.rep <- NA
		shat.pure.error <- NA
		out$pure.ss <- 0
		YM <- y
		weightsM <- wt
		xM <- x[!duplicated(rep.info)]
	}
	else {
		##
		## do a simple 1-way ANOVA to get the replication error
		##
		rep.info.aov <- fast.1way(rep.info, y, wt)
		shat.pure.error <- sqrt(rep.info.aov$MSE)
		shat.rep <- shat.pure.error
		YM <- rep.info.aov$means
		weightsM <- rep.info.aov$w.means
		xM <- as.matrix(x[!duplicated(rep.info)])
		out$pure.ss <- rep.info.aov$SSE
		if(verbose) {
			cat(" rep info", fill = TRUE)
			print(rep.info.aov)
		}
	}
	#
	# over write the real x's and y's with the collapsed group means. 
	# these orignal values are called xraw yraw. 
	#
	out$np <- length(xM)
	out$y <- c(YM)
	out$x <- c(xM)
	out$wt <- weightsM
	out$shat.rep <- shat.rep
	out$shat.pure.error <- shat.pure.error
	xgrid <- sort(xM)
	out$trace <- NA
	#
	# find lambda's if df's are given
	#
	if(!is.na(df[1])) {
		lam <- rep(0, length(df))
		for(k in 1:length(df)) {
			lam[k] <- sreg.df.to.lambda(df[k], xM, weightsM)
		}
	}
	if(do.cv) {
		a <- gcv.sreg(out, lambda = lam, cost = cost, offset = offset,
			nstep.cv = nstep.cv, verbose = verbose, trmin = trmin,
			trmax = trmax, find.min = find.min, method = method,
			rmse = rmse)
		# if the spline should be evaluated at the GCV solutionwipe out lam grid
		# and just use GCV lambda.
		out$gcv.grid <- a$gcv.grid
		out$lambda.est <- a$lambda.est
		# 
		# save  GCV estimate is that is what is needed
		if(eval.cv) {
			lam <- a$lambda.est[method, "lambda"]
			out$trace <- a$lambda.est[method, "trA"]
			out$shat.GCV <- a$lambda.est[method, "shat"]
		}
		#
		# save trace and estimate of sigma an the traces
		if(!find.min) {
			out$trace <- c(a$gcv.grid[, "trA"])
			out$shat.GCV <- c(a$gcv.grid[1, "shat"])
		}
	}
	b <- list()
	# lam can either be  a grid or just the GCV value 
	NL <- length(lam)
	NG <- length(xgrid)
	h <- log(lam)
	residuals <- matrix(NA, ncol = NL, nrow = length(out$yraw))
	job <- as.integer(c(0, 3, 0))
	predicted <- matrix(NA, ncol = NL, nrow = NG)
	if(find.diagA) {
		diagA <- matrix(NA, ncol = NL, nrow = out$np)
		# add switch to find diag of A. 
		job <- as.integer(c(3, 3, 0))
	}
	for(k in 1:NL) {
		#
		# call cubic spline, note lambda is passed in log scale. 
		# spline solution evaluated at xgrid
		# 
		b <- .Fortran("css",
			h = as.double(h[k]),
			npoint = as.integer(out$np),
			x = as.double(out$x),
			y = as.double(out$y),
			wt = as.double(1/sqrt(out$wt)),
			sy = as.double(rep(0, out$np)),
			trace = as.double(0),
			diag = as.double(c(cost, offset, rep(0, (out$np - 2)))),

				cv = as.double(0),
			ngrid = as.integer(NG),
			xg = as.double(xgrid),
			yg = as.double(rep(0, NG)),
			job = as.integer(job),
			ideriv = as.integer(0),
			ierr = as.integer(0), PACKAGE="fields")
		#
		if(verbose) {
			print(c(b$lambda, b$trace))
		}
		if(find.diagA) {
			diagA[, k] <- b$diag
		}
		residuals[, k] <- out$yraw - splint(out$x, b$sy, out$xraw)
		predicted[, k] <- b$yg
	}
	out$call <- call
	out$lambda <- lam
	out$do.cv <- do.cv
	out$residuals <- residuals
	out$fitted.values <- out$yraw - residuals
	out$predicted <- list(x = xgrid, y = predicted)
	if(length(lambda[1]) == 1) {
		out$eff.df <- out$trace[1]
	}
	if(find.diagA) {
		out$diagA <- diagA
	}
	out
}
