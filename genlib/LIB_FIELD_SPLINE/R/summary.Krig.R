"summary.Krig" <-
function (object, digits = 4, ...) 
{
    x <- object

# lambda est may not be available if lambda has been supplied by user.

if (!is.na(x$lambda.est[1])){
         l.est<- x$lambda.est}
else{
         l.est<- NA} 

    summary <- list(call = x$call, num.observation = length(x$residuals), 
        enp = x$eff.df, nt = x$nt, df.drift= sum(x$ind.drift),
             res.quantile = quantile(x$residuals, 
            seq(0, 1, 0.25)), shat.MLE = x$shat.MLE, shat.GCV = x$shat.GCV, 
        rhohat = x$rhohat, m = x$m, lambda = x$lambda, cost = x$cost, 
        rho = x$rho, sigma2 = x$sigma2, 
        num.uniq = length(x$yM), knot.model = x$knot.model, np = x$np, 
        method = x$method, lambda.est = l.est, 
        shat.pure.error = x$shat.pure.error, args=x$args)
    class(summary) <- "summary.Krig"
    summary$covariance <- cor(x$fitted.values * sqrt(x$weights), 
        (x$y) * sqrt(x$weights))^2
    hold <- (sum((x$y - mean(x$y))^2) - sum(x$residuals^2))/(sum((x$y - 
        mean(x$y))^2))
    summary$adjr2 <- 1 - ((length(x$residuals) - 1)/(length(x$residuals) - 
        x$eff.df)) * (1 - hold)
    summary$digits <- digits
    summary$cov.function <- as.character(x$cov.function.name)
    summary$correlation.model <- x$correlation.model
    summary$sum.gcv.lambda <- summary.gcv.Krig(x, x$lambda)
    summary
}

