"Krig.parameters" <-
function (obj, mle.calc = obj$mle.calc) 
{


# if nondiag W is supplied then use it.
# otherwise assume a diagonal set of weights.
#
# NOTE: calculation of  shat involves full set of obs
# not those colllapsed to the mean.
  
  if(obj$nondiag.W){
    shat.GCV <- sqrt( 
          sum( (obj$W2%d*% obj$residuals)^2)/(length(obj$y) - obj$eff.df)
              )}
  else{
          shat.GCV <- sqrt(
             sum(
                 (obj$weights*obj$residuals^2)/
                                  (length(obj$y) - obj$eff.df)
                        ))}

  
    if (mle.calc) {
        rhohat <- sum(obj$c * obj$yM)/(obj$N - obj$nt)
        shat.MLE <- sqrt(rhohat * obj$lambda)
    }
    else {
        rhohat <- shat.MLE <- NA
    }
    list(shat.GCV = shat.GCV, shat.MLE = shat.MLE, rhohat = rhohat)
}

