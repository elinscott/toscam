\name{predict.se}
\alias{predict.se}
\title{
  Standard errors of predictions 
}
\description{
Calculates the standard error of predictions. This is
usually the fitted object from a function estimate such as from Krig
or Tps. 
}
\usage{
predict.se(object, ...)
}
\arguments{
\item{object}{
A fitted model object of a certain class 
}
\item{\dots}{
Additional arguments to be passed to 
a particular method. e.g. a grid.list or model specification. }
}
\value{
A vector of standard errors for the predicted values. 
}
\details{
This function is generic and will call the appropriate function to calculate 
the standard errors for the object class. The prediction standard error is
for
the estimated function or parameters (a mean value) not for the
prediction of a new observation.
}
\seealso{
predict, predict.surface.se, predict.se.Krig  
}
\examples{
  fit<-Krig(ozone$x,ozone$y, Covariance="Matern",theta=50, smoothness=1)
  predict.se(fit)                             # std errors of predictions 
#
# create a grid of points  
  xg<- make.surface.grid(
        list(East.West=seq(-15,15,,20),North.South=seq(-20,20,,20) ) ) 
  out<- predict.se(fit,xg)
  image.plot( as.surface( xg, out))

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
