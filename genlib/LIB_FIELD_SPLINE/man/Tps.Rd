\name{Tps}
\alias{Tps}
\title{
  Thin plate spline regression 
}
\description{
Fits a thin plate spline surface to irregularly spaced data. The 
smoothing parameter is chosen by generalized cross-validation. The assumed 
model is additive  Y = f(X) +e  where f(X) is a d dimensional surface. 
This is a special case of the spatial process estimate.   
}
\usage{
Tps(x, Y, m = NULL, p = NULL, scale.type = "range", ...)
}
\arguments{
To be helpful, a more complete list of arguments are described that are the 
same as those for the Krig function.  

\item{x}{
Matrix of independent variables. Each row is a location or a set of 
independent covariates. 
}
\item{Y}{
Vector of dependent variables. 
}
\item{m}{
A polynomial function of degree (m-1) will be  
included in the model as the drift (or spatial trend) component. 
Default is the value such that 2m-d is greater than zero where d is the 
dimension of x. 
}

\item{p}{
Exponent for radial basis functions. Default is 2m-d. 
}

\item{scale.type}{
The independent variables and knots are scaled to the specified 
scale.type. 
By default the scale type is "range", whereby 
the locations are transformed  
 to the interval (0,1) by forming (x-min(x))/range(x) for each x. 
Scale type of "user" allows specification of an x.center and x.scale by 
the 
user. The default for "user" is mean 0 and standard deviation 1. Scale 
type of "unscaled" does not scale the data.  
}

\item{\dots}{
Any argument that is valid for the Krig function. Some of the main ones
are listed below. 

\describe{

\item{lambda}{
Smoothing parameter that is the ratio of the error variance (sigma**2) 
to the scale parameter of the  
covariance function. If omitted this is estimated by GCV. 
}

\item{df}{
The effective number of parameters for the fitted surface. Conversely,
N- df, where N is the total number of observations is the degrees of
freedom associated with the residuals.
This is an alternative to specifying lambda and much more interpretable.
}

\item{cost}{
Cost value used in GCV criterion. Corresponds to a penalty for
increased number of parameters. The default is 1.0 and corresponds to the
usual GCV.
}


\item{weights}{
Weights are proportional to the reciprocal variance of the measurement  
error. The default is no weighting i.e. vector of unit weights. 
}
\item{nstep.cv}{
Number of grid points for minimum GCV search. 
}
\item{x.center}{
Centering values are subtracted from each column of the x matrix. 
Must 
have scale.type="user".
}
\item{x.scale}{
Scale values that divided into each column after centering. 
Must 
have scale.type="user".
}
\item{rho}{
Scale factor for covariance. 
}
\item{sigma2}{
Variance of errors or if weights are not equal to 1 the variance is
sigma**2/weight. 
}

\item{method}{
Determines what "smoothing" parameter should be used. The default
is to estimate standard GCV
Other choices are: GCV.model, GCV.one, RMSE, pure error and REML. The
differences are explained below.
}
\item{verbose}{
If true will print out all kinds of intermediate stuff.  
}
\item{mean.obj}{
Object to predict the mean of the spatial process. 
}
\item{sd.obj}{
Object to predict the marginal standard deviation of the spatial process. 
}

\item{null.function}{
An R function that creates the matrices for the null space model.
The default is fields.mkpoly, an R function that creates a polynomial
regression matrix with all terms up to degree m-1. (See Details)
}

\item{offset}{
The offset to be used in the GCV criterion. Default is 0. This would be 
used when Krig/Tps is part of a backfitting algorithm and the offset has
to be included  to reflect other model degrees of freedom. 
}
}
}
}
\value{
A list of class Krig. This includes the predicted surface of 
fitted.values and the residuals. The results of the grid 
search minimizing the generalized cross validation function is 
returned in gcv.grid. 
Please see the documentation on Krig for details of the returned 
arguments.  
}
\details{
A thin plate spline is result of minimizing the residual sum of 
squares subject to a constraint that the function have a certain 
level of smoothness (or roughness penalty). Roughness is 
quantified by the integral of squared m-th order derivatives. For one 
dimension and m=2 the roughness penalty is the integrated square of 
the second derivative of the function. For two dimensions the 
roughness penalty is the integral of  

        (Dxx(f))**22 + 2(Dxy(f))**2 + (Dyy(f))**22 

(where Duv denotes the second partial derivative with respect to u 
and v.) Besides controlling the order of the derivatives, the value of 
m also determines the base polynomial that is fit to the data. 
The degree of this polynomial is (m-1). 

The smoothing parameter controls the amount that the data is 
smoothed. In the usual form this is denoted by lambda, the Lagrange 
multiplier of the minimization problem. Although this is an awkward 
scale, lambda =0 corresponds to no smoothness constraints and the data 
is interpolated.  lambda=infinity corresponds to just fitting the 
polynomial base model by ordinary least squares.  

This estimator is implemented by passing the right generalized covariance
function based on radial basis functions to the more general function
Krig.  One advantage of this implementation is that once a Tps/Krig object
is created the estimator can be found rapidly for other data and smoothing
parameters provided the locations remain unchanged. This makes simulation
within R efficient (see example below). Tps does not currenty support the
knots argument where one can use a reduced set of basis functions. This is
mainly to simplify and a good alternative using knots would be to use a
valid covariance from the Matern family and a large range parameter. 

 }

\section{References}{
See "Nonparametric Regression and Generalized Linear Models"  
by Green and Silverman. 
See "Additive Models" by Hastie and Tibshirani. 
}
\seealso{
Krig, summary.Krig, predict.Krig, predict.se.Krig, plot.Krig,  
\code{\link{surface.Krig}}, 
\code{\link{sreg}} 
}
\examples{
#2-d example 

fit<- Tps(ozone$x, ozone$y)  # fits a surface to ozone measurements. 

set.panel(2,2)
plot(fit) # four diagnostic plots of  fit and residuals. 
set.panel()

summary(fit)

# predict onto a grid that matches the ranges of the data.  

out.p<-predict.surface( fit)
image( out.p) 
surface(out.p) # perspective and contour plots of GCV spline fit 
# predict at different effective 
# number of parameters 
out.p<-predict.surface( fit,df=10)

#1-d example 
out<-Tps( rat.diet$t, rat.diet$trt) # lambda found by GCV 
plot( out$x, out$y)
lines( out$x, out$fitted.values)

# 
# compare to the ( much faster) one spline algorithm 
#  sreg(rat.diet$t, rat.diet$trt) 
# 

# Adding a covariate to the fixed part of model
# Note: this is a fairly big problem numerically (850+ locations)

Tps( RMprecip$x,RMprecip$y, Z=RMprecip$elev)-> out

surface( out, drop.Z=TRUE)
US( add=TRUE, col="grey")

#
# simulation reusing Tps/Krig object
#
fit<- Tps( rat.diet$t, rat.diet$trt)
true<- fit$fitted.values
N<-  length( fit$y)
temp<- matrix(  NA, ncol=50, nrow=N)
sigma<- fit$shat.GCV
for (  k in 1:50){
ysim<- true + sigma* rnorm(N) 
temp[,k]<- predict(fit, y= ysim)
}
matplot( fit$x, temp, type="l")


# 
#4-d example 
fit<- Tps(BD[,1:4],BD$lnya,scale.type="range") 

# plots fitted surface and contours 
# default is to hold 3rd and 4th fixed at median values 

surface(fit)   



}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
