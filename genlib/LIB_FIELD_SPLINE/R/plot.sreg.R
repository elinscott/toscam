"plot.sreg" <-
function(x, digits = 4, which = 1:4,...)
{
        out<- x # hack S3

	if(any(which==1)) {
		plot(out$xraw, out$yraw, ylab = "predicted", xlab = " X", bty
			 = "n", ...)
		matlines(out$predicted$x, out$predicted$y, lty = 1)
	}
	if(any(which==2) & length(out$lambda) == 1) {
		plot(out$fitted.values, out$residuals, ylab = "residuals",
			xlab = " predicted values", bty = "n", ...)
		yline(0)
	}
	if(any(which==3)) {
		if(nrow(out$gcv.grid) > 1) {
			# trim off + infinity due to pole in the denominator of GCV function
			#with cost
			ind <- out$gcv.grid[, 3] < 1e+19
			out$gcv.grid <- out$gcv.grid[ind,  ]
			yr <- range(unlist(out$gcv.grid[, 3:5]), na.rm=TRUE)
	plot(out$gcv.grid[, 2], out$gcv.grid[, 3], xlab = 
				"Eff. parameters", ylab = " GCV function",
				bty = "n", ylim = yr, log = "y", ...)
			lines(out$gcv.grid[, 2], out$gcv.grid[, 4], lty = 2)
			lines(out$gcv.grid[, 2], out$gcv.grid[, 5], lty = 1)
			xline(out$eff.df)
			title("GCV-points , solid- GCV model,\ndashed- GCV one",
				cex = 0.6)
		}
	}
	if(any(which==4)) {
		if(length(out$lambda) == 1) {
			hist(out$residuals, xlab = "Residuals", main="")
		}
		else {
			bplot(out$residuals, labels=
                                    format(round(out$trace,1)), 
				xlab = "eff df", srt=90)

		title("Residuals")
		}
	}
}
