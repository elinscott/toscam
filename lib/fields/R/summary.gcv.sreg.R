"summary.gcv.sreg" <-
function(object, lambda, cost = 1, nstep.cv = 20, offset = 0, verbose = 
TRUE,...)
{
        out<- object # hack S3
	shat.pure.error <- out$shat.pure.error
	pure.ss <- out$pure.ss
	nt <- 2
	np <- out$np
	N <- out$N
	out$cost <- cost
	out$offset <- offset
	lambda.est <- rep(NA, 6)
	names(lambda.est) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model",
		"shat")
	#
	# fill in stuff for this  lambda
	lambda.est[1] <- lambda
	temp <- sreg.fit(lambda, out)
	lambda.est[2] <- temp$trace
	lambda.est[3] <- temp$gcv
	lambda.est[4] <- temp$gcv.one
	if(!is.na(shat.pure.error)) {
		lambda.est[5] <- temp$gcv.model
	}
	lambda.est[6] <- temp$shat
	lambda.est
}
