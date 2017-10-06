"Krig.updateY" <-
function(out, Y, verbose = FALSE, yM=NA )
{
	#given new Y values but keeping everything else the same finds the 
	#new u vector and pure error SS associated with the Kriging estimate
	# the steps are
	# 1) standardize if neccesary
	# 2) find means, in the case of replicates
	# 3) based on the decomposition, multiply a weighted version of yM
	#    with a large matrix extracted from teh Krig object out. 
	#
	# The out object may be large. This function is written so that out is # #not changed with the hope that it is not copied locally  in this  #function .
	# All of the output is accumulated in the list out2
	#STEP 1
	#
	# transform Y by mean and sd if needed
	#
	if(out$correlation.model) {
		Y <- (Y - predict(out$mean.obj, out$x))/predict(out$sd.obj,
			out$x)
		if(verbose)
			print(Y)
	}
	#
	#STEP 2
        if( is.na( yM[1])){
          out2 <- Krig.ynew(out, Y)}
        else{
          out2<- list( yM= yM, 
          shat.rep = NA,shat.pure.error = NA, pure.ss = NA)}

	if(verbose) {
		print(out2)
	}
	#
	#STEP3
	#
	# Note how matrices are grabbed from the Krig object
	#
	if(verbose) cat("Type of decomposition", out$decomp, fill = TRUE)

	if(out$decomp == "DR") {
		#
		#
		u <- t(out$matrices$G) %*% t(out$matrices$X) %*% (out$weightsM *
			out2$yM)
		#
		# find the pure error sums of sqaures.  
		#
		temp <- out$matrices$X %*% out$matrices$G %*% u
		temp <- sum((out$W2%d*%(out2$yM - temp) )^2)
		out2$pure.ss <- temp + out2$pure.ss
		if(verbose) {
			cat("pure.ss", fill = TRUE)
			print(temp)
			print(out2$pure.ss)
		}
	}
	#####
	##### end DR decomposition block 
	#####
	####
	#### begin WBW decomposition block
	####
	if(out$decomp == "WBW") {
		#### decomposition of Q2TKQ2
		u <- c(rep(0, out$nt), t(out$matrices$V) %*% qr.q2ty(out$
			matrices$qr.T, out$W2%d*%out2$yM ))
		if(verbose)
			cat("u", u, fill = TRUE)
		#
		# pure error in this case from 1way ANOVA 
		#
		if(verbose) {
			cat("pure.ss", fill = TRUE)
			print(out2$pure.ss)
		}
	}
	#####
	##### end WBW block
	#####
	out2$u <- u
	out2
}
