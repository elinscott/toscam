"predict.sreg" <-
function(object, x, derivative = 0, model = 1,...)
{
        
	if(missing(x))
		x <- object$xraw
	c(
    splint(object$predicted$x, object$predicted$y[, model], x, 
             derivative = derivative)
)
}
