\name{splint}
\alias{splint}
\title{
  Cubic spline interpolation 
}
\description{
A fast, FORTRAN based function for cubic spline interpolation. 
}
\usage{
splint(x, y, xgrid, derivative=0)
}
\arguments{
\item{x}{
The x values that define the curve or a two column matrix of  
x and y values. 
}
\item{y}{
The y values that are paired with the x's. 
}
\item{xgrid}{
The grid to evaluate the fitted cubic interpolating curve. 
}
\item{derivative}{
Indicates whether the function or a a first or second derivative 
should be evaluated. 
}
}
\value{
A vector consisting of the spline evaluated at the grid values. 
}
\details{
Fits a piecewise interpolating cubic polynomial to the x and y values
assuming that the second and third derivatives are zero at the range of the x 
values. For this reason extrapolation outside the range of the x 
values will be a linear function. 
}
\section{References}{
See Additive Models by Hastie and Tibshriani. 
}
\seealso{
sreg,  Tps  
}
\examples{
x<- seq( 0, 120,,200)

splint(rat.diet$t, rat.diet$trt,x )-> y

plot( rat.diet$t, rat.diet$trt)
lines( x,y)
#( this is weird!)
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
