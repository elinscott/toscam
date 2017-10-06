"plot.Krig" <-
function (x, digits = 4, which = 1:4,  ...) 
{
    out <- x
#
#   don't do plots 2:4 if a fixed lambda
#
    if( x$fixed.model){ which <- 1 }

    fitted.values <- predict(out)
    std.residuals <- (out$residuals * sqrt(out$weights))/out$shat.GCV



    if ( any(which==1) ) {
        temp <- summary(out)
        plot(fitted.values, out$y, ylab = "Y", xlab = " predicted values", 
            bty = "n", ...)
        abline(0, 1)
        hold <- par("usr")
        text(hold[1], hold[4], paste(" R**2 = ", format(round(100 * 
            temp$covariance, 2)), "%", sep = ""), cex = 0.8, 
            adj = 0)
    }
    if (any(which==2) )  {
        plot(fitted.values, std.residuals, ylab = "(STD) residuals", 
            xlab = " predicted values", bty = "n", ...)
        yline(0)
        hold <- par("usr")
        text(hold[1], hold[4], 
        paste(" RMSE =", 
            format(signif(sqrt(sum(out$residuals^2)/(temp$num.observation - temp$enp)), digits))), cex = 0.8, adj = 0)
    }
    if (any(which==3) ) {
        if (nrow(out$gcv.grid) > 1) {
            ind <- out$gcv.grid[, 3] < 1e+19
            out$gcv.grid <- out$gcv.grid[ind, ]
            yr <- range(unlist(out$gcv.grid[, 3:5]), na.rm = TRUE)
            plot(out$gcv.grid[, 2], out$gcv.grid[, 3], xlab = "Eff. number of parameters", 
                ylab = " GCV function", bty = "n", ylim = yr, 
                log = "y", ...)

            lines(out$gcv.grid[, 2], out$gcv.grid[, 4], lty = 2)
            lines(out$gcv.grid[, 2], out$gcv.grid[, 5], lty = 1)
            xline(out$eff.df)
            title("GCV-points , solid- GCV model,\ndashed- GCV one", 
                cex = 0.6)
    box()
        }
    }

    if (any(which==4)) {
        hist(std.residuals)
    }

}

