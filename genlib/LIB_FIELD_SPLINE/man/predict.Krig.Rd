\name{predict.Krig}
\alias{predict.Krig}
\title{
  Evaluation of Krig spatial process estimate.  
}
\description{
Provides predictions from the Krig spatial process estimate at arbitrary
points, new data (Y) or other values of the smoothing parameter (lambda)
including a GCV estimate. 
}
\usage{
predict.Krig(
object, x = NULL, Z = NULL, drop.Z = FALSE, just.fixed
                 = FALSE, lambda = NA, df = NA, model = NA,
                 eval.correlation.model = TRUE, y = NULL, yM = NULL,
                 verbose = FALSE, ...)

}
\arguments{
\item{object}{
Fit object from the Krig or Tps function. 
}
\item{x}{
Matrix of x values on which to evaluate the kriging surface. 
If omitted, the data x values, i.e. out\$x  will be used. 
}
\item{Z}{ Vector/Matrix of additional covariates to be included in 
fixed part of spatial model}
\item{drop.Z}{ If TRUE only spatial fixed part of model is evaluated. 
i.e. Z covariates are not used. 
}
\item{just.fixed}{
Only fixed part of model is evaluated}

\item{lambda}{
Smoothing parameter. If omitted, out\$lambda will be used. 
(See also df and gcv arguments)
}
\item{df}{
Effective degrees of freedom for the predicted surface. This can be used
in place of lambda ( see the function Krig.df.to.lambda)
}
\item{model}{
Generic argument that may be used to pass a different lambda.  
}
\item{eval.correlation.model}{
If true ( the default) will multiply the predicted function by marginal
sd's
and add the mean function. This usually what one wants. If false will
return predicted surface in the standardized scale. The main use of this
option is a call from Krig to find MLE's of rho and sigma2     
}
\item{y}{
Evaluate the estimate using the new data vector y (in the same order as
the old data). This is equivalent
to recomputing the Krig object with this new data but is more efficient
because many pieces can be reused. Note that the x values are assumed to
be the same. 
}
\item{yM}{
If not NULL evaluate the estimate using this vector as the replicate mean 
data. That is, assume the full data has been collapsed into replicate 
means in the same order as xM. The replicate weights are assumed to be 
the same as the original data. (weightsM)
}

 \item{verbose}{
Print out all kinds of intermediate stuff for debugging 
}

\item{\dots}{Other arguments passed to predict.}
}
\value{
Vector of predicted responses
}
\details{
The main goal in this function is to reuse the Krig object to rapidly
evaluate different estimates. Thus there is flexibility in changing the
value of lambda and also the independent data without having to recompute
the matrices associated with the Krig object. The reason this is possible
is that most on the calculations depend on the observed locations not on
lambda or the observed data.  
}
\seealso{
Krig, predict.surface gcv.Krig 
}

\examples{
Krig(ozone$x,ozone$y, theta=50) ->fit
predict( fit) # gives predicted values at data points

# only the fixed part of the model

predict( fit, just.fixed=TRUE) 

# in this case the default is a linear spatial drift  (m=2) and there 
# are no additional covariates

grid<- make.surface.grid( list( seq( -40,40,,15), seq( -40,40,,15)))

look<- predict(fit,grid) # evaluate on a grid of points

# some useful graphing functions
out.p<- as.surface( grid, look) # reformat into $x $y $z image-type object
contour( out.p)  


# refit with 10 degrees of freedom in surface

look<- predict(fit,grid, df=15)

# refit with random data 

look<- predict( fit, grid, y= rnorm( 20))

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
