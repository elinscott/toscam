\name{Krig}
\alias{Krig}
\alias{resid.Krig}
\alias{fitted.Krig}
\alias{coef.Krig}

\title{
  Kriging surface estimate  
}
\description{
Fits a surface to irregularly spaced data. The Kriging model assumes  
that the unknown  function is a realization  of a Gaussian 
random spatial processes. The assumed model is additive 
Y = P(x) +  Z(X) + e, where P is a low order polynomial and Z is a 
mean zero, 
Gaussian stochastic process with a  
covariance that is unknown up to a scale constant. The main advantages of
this function are the flexibility in specifying the covariance as an
R language function and also the supporting functions plot, predict,
predict.se, surface for
subsequent analysis. Krig also supports a correlation model where the mean
and marginal variances are supplied. 
}
\usage{
Krig(
x, Y, cov.function = "stationary.cov", lambda = NA, df
                 = NA, Z = NULL, cost = 1, knots = NA, weights = NULL,
                 m = 2, nstep.cv = 80, scale.type = "user", x.center =
                 rep(0, ncol(x)), x.scale = rep(1, ncol(x)), rho = NA,
                 sigma2 = NA, method = "GCV", verbose = FALSE, mean.obj
                 = NA, sd.obj = NA, null.function =
                 "Krig.null.function", wght.function = NULL, offset =
                 0, outputcall = NULL, na.rm = TRUE, cov.args = NULL,
                 chol.args = NULL, null.args = NULL, wght.args = NULL,
                 W = NULL, give.warnings = TRUE, ...)

fitted.Krig(object,...)

coef.Krig(object,...)

resid.Krig(object,...)

}
\arguments{
\item{x}{
Matrix of independent variables. These could the locations for spatial
data or the indepedent variables in a regression. 
}
\item{Y}{
Vector of dependent variables. These are the values of the surface
(perhaps with measurement error) at the locations or the dependent
response in a regression.
}
\item{cov.function}{
Covariance function for data in the form of an R function (see 
Exp.simple.cov as an example).  
Default assumes that correlation is an exponential function of distance.
See also \code{stationary.cov} for more general choice of covariance 
shapes. \code{exponential.cov} will be faster if only the exponential 
covariance form is needed. 
}
\item{Z}{
A vector of matrix of covariates to be include in the fixed part of the 
model. If NULL (default) no addtional covariates are included.}

\item{lambda}{
Smoothing parameter that is the ratio of the error variance (sigma**2) 
to the scale parameter of the  
covariance function (rho). If omitted this is estimated by GCV ( see
method below). 
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

\item{knots}{
A matrix of locations similar to x. These can define an alternative set of
basis functions for representing the estimate. One choice may be a
space-filling subset of the original x locations, thinning out the
design where locations cluster. The 
default is to put a "knot" at all unique locations. (See details.)
}

\item{weights}{
Weights are proportional to the reciprocal variance of the measurement  
error. The default is equal weighting i.e. vector of unit weights. 
}

\item{m}{
A polynomial function of degree (m-1) will be  
included in the model as the drift (or spatial trend) component. 
The "m" notation is from thin-plate splines where m is the 
derivative in the penalty function. With m=2 as the default a linear 
model in the locations  will be fit a fixed part of the model. 
}

\item{nstep.cv}{
Number of grid points for the coarse grid search  to 
minimize the GCV  RMLE and other related criterian for finding lambda. 
}

\item{scale.type}{
This is a character string among: "range", "unit.sd", "user", "unscaled".
The independent variables and knots are scaled to the specified scale.type. 
By default no scaling is done. This usuall makes sense for spatial locations.
Scale type of "range" scales the data to the interval (0,1) by forming 
(x-min(x))/range(x) for each x. Scale type of "unit.sd" 
Scale type of "user" allows specification of an x.center and x.scale by the 
user. The default for "user" is mean 0 and standard deviation 1. Scale 
type of "unscaled" does not scale the data. 
}

\item{x.center}{
Centering values to be subtracted from each column of the x matrix. 
}

\item{x.scale}{
Scale values that are divided into each column after centering. 
}

\item{rho}{
Scale factor for covariance. 
}

\item{sigma2}{
Variance of the errors, often called the nugget variance. If weights are
specified then the error variance is sigma2 divided by weights. 
Note that lambda is defined as the ratio  sigma2/rho.
}

\item{method}{
Determines what "smoothing" parameter should be used. The default 
is to estimate standard GCV 
Other choices are: GCV.model, GCV.one, RMSE, pure error and REML. The 
differences are explained below.  
}

\item{verbose}{
If true will print out all kinds of intermediate stuff. Default is false,
of course as this is used mainly for debugging.
}

\item{mean.obj}{
Object to predict the mean of the spatial process. This used in when
fitting a correlation model with varying spatial means and varying
marginal variances. (See details.)
}

\item{sd.obj}{
Object to predict the marginal standard deviation of the spatial process. 
}

\item{null.function}{
An R function that creates the matrices for the null space model.  
The default is fields.mkpoly, an R function that creates a polynomial 
regression matrix with all terms up to degree m-1. (See Details)  
}

\item{wght.function}{
An R function that creates a weights matrix to the observations.  
This is only needed if the weight matirx has off diagonal elements. 
The default is NULL indicating that the weight matrix is a diagonal, based
on the weights argument. (See details)
}

\item{offset}{
The offset to be used in the GCV criterion. Default is 0. This would be 
used when Krig is part of a backfitting algorithm and the offset is 
other model degrees of freedom from other regression components. 
}

\item{outputcall}{
If NULL the output object will have a \$call argument based on this call. 
If no NULL the output call will have whatever is passed. This is kludge
for the Tps function so that it return a Krig object but have the right
call argument. Sorry no one promised that fields would be pretty. 
}

\item{cov.args}{
A list with the arguments to call the covariance function. (in addition to the locations)
}

\item{na.rm}{If TRUE NAs will be removed from the \code{y} vector and the 
corresponding rows of \code{x} -- with a warning. 
If FALSE Krig will just stop with a message. Once NAs are removed all 
subsequent analysis in fields does not use those data. }

\item{chol.args}{
Arguments to be passed to the cholesky decomposition in Krig.engine.fixed. 
The default if NULL, assigned at the top level of this function, is 
list( pivot=FALSE). This argument is useful when working with 
the sparse matrix package. }

\item{wght.args}{
Optional arguments to be passed to the weight function  (wght.function) 
used to create the observation weight matrix.} 

\item{W}{The explicit observatoin weight matrix.}

\item{null.args}{ Extra arguments for the null space function 
\code{null.function}. 
If \code{fields.mkpoly} is passed as 
\code{null.function} then this is set to a list with the
value of \code{m}. So the default is use a polynomial of degree m-1 for the
null space (fixed part) of the model. 
}

\item{give.warnings}{ If TRUE warnings are given in gcv grid search limits.
If FALSE warnings are not given. Best to leave this TRUE!}

\item{\dots}{
 Optional arguments that appear are assumed to be additional arguments
to the covariance function. Or are included in methods functions (resid,
fitted, coef) as a
required argument.}

\item{object}{ A Krig object}



}
\value{
A object of class Krig. This includes the predicted values in  
fitted.values and the residuals in residuals. The results of the grid 
search to minimize the generalized cross validation function are 
returned in gcv.grid. 

The coef.Krig function only returns the coefficients, "d", associated with the 
fixed part of the model (also known as the null space or spatial drift).

\item{call}{
Call to the function 
}
\item{y}{
Vector of dependent variables. 
}
\item{x}{
Matrix of independent variables. 
}
\item{weights}{
Vector of weights. 
}
\item{knots}{
Locations used to define the basis functions.  
}
\item{transform}{
List of components used in centering and scaling data. 
}
\item{np}{
Total number of parameters in the model. 
}
\item{nt}{
Number of parameters in the null space. 
}
\item{matrices}{
List of matrices from the decompositions (D, G, u, X, qr.T). 
}
\item{gcv.grid}{
Matrix of values from the GCV grid search. The first column 
is the grid of lambda values used in the search, the second column  
is the trace of the A matrix, the third column is the GCV values and 
the fourth column is the estimated value of sigma conditional on the vlaue
of lambda.  
}
\item{lambda.est}{
A table of estimated smoothing parameters with corresponding degrees 
of freedom and estimates of sigma found by different methods.  
}
\item{cost}{
Cost value used in GCV criterion. 
}
\item{m}{
Order of the polynomial space: highest degree polynomial is (m-1). 
This is a fixed part of the surface often referred to as the drift 
or spatial trend.  
}
\item{eff.df}{
Effective degrees of freedom of the model. 
}
\item{fitted.values}{
Predicted values from the fit. 
}
\item{residuals}{
Residuals from the fit. 
}
\item{lambda}{
Value of the smoothing parameter used in the fit. 
}
\item{yname}{
Name of the response. 
}
\item{cov.function}{
Covariance function of the model. 
}
\item{beta}{
Estimated coefficients in the ridge regression format 
}
\item{d}{
Estimated coefficients for the polynomial basis functions that span the 
null space 
}
\item{fitted.values.null}{
Fitted values for just the polynomial part of the estimate 
}
\item{trace}{
Effective number of parameters in model. 
}
\item{c}{
Estimated coefficients for the basis functions derived from the 
covariance. 
}
\item{coefficients}{
Same as the beta vector. 
}
\item{just.solve}{
Logical describing if the data has been interpolated using the basis  
functions.  
}
\item{shat}{
Estimated standard deviation of the measurement error (nugget effect). 
}
\item{sigma2}{
Estimated variance of the measurement error (shat**2). 
}
\item{rho}{
Scale factor for covariance.  COV(h(x),h(x\code{)) = rho*cov.function(x,x}) 
If the covariance is actually a 
correlation function then rho is also the "sill". 
}
\item{mean.var}{
Normalization of the covariance function used to find rho. 
}
\item{best.model}{
Vector containing the value of lambda, the estimated variance of the  
measurement error and the scale factor for covariance used in the fit. 
}
}
\details{

This function produces a object of class Krig. With this object it is
easy to subsequently predict with this fitted surface, find standard
errors, alter the y data ( but not x), etc. 

The Kriging model is: Y.k= P(x.k) + Z(x.k) + e.k

where ".k" means subscripted by k, Y is the dependent variable observed
at location x.k, P is a low order polynomial, Z is a mean zero, Gaussian
field with covariance function K and e is assumed to be independent
normal errors. The estimated surface is the best linear unbiased
estimate (BLUE) of P(x) + Z(x) given the observed data. For this
estimate K, is taken to be rho*cov.function and the errors have variance
sigma**2. In more conventional geostatistical terms rho is the "sill" if
the covariance function is actually a correlation function and sigma**2
is the nugget variance or measure error variance (the two are confounded
in this model.)  If the weights are given then the variance of e.k is
sigma**2/ weights.k . In the case that the weights are specified as a
matrix, W,  using the wght.function option then the assumed covariance matrix
for the errors is sigma**2 Wi, where Wi is the inverse of W. 

If these parameters rho and sigma2 are omitted in the call, then they are
estimated in the following way. If lambda is given, then sigma2 is
estimated from the residual sum of squares divided by the degrees of
freedom associated with the residuals.  Rho is found as the difference
between the sums of squares of the predicted values having subtracted off
the polynomial part and sigma2.  

A useful extension of a stationary correlation to a nonstationary
covariance is what we term a correlation model. 
If mean and marginal standard deviation objects are included in the call.  
Then the observed data is standardized based on these functions.  The
spatial process is then estimated with respect to the standardized scale.
However for predictions and standard errors the mean and standard
deviation surfaces are used to produce results in the original scale of
the observations.

The GCV function has several alternative definitions when replicate
observations are present or if one uses a reduced set knots.  Here are the
choices based on the method argument:  

GCV: leave-one-out GCV. But if
there are replicates it is leave one group out. (Wendy and Doug prefer
this one.)  

GCV.one: Really leave-one-out GCV even if there are replicate
points.  This what the old tps function used in FUNFITS.
 
rmse: Match the estimate of sigma**2 to a external value ( called rmse)  

pure error: Match the estimate of sigma**2 to the estimate based on
replicated data (pure error estimate in ANOVA language).  

GCV.model:
Only considers the residual sums of squares explained by the basis
functions.  

WARNING: The covariance functions often have a nonlinear parameter(s) that
often control the strength of the correlations as a function of separation,
usually referred to as the range parameter. This parameter must be
specified in the call to Krig and will not be estimated. 

}
\section{References}{
See "Additive Models" by Hastie and Tibshirani, "Spatial Statistics" by    
Cressie and the FIELDS manual. 
}
\seealso{
summary.Krig, predict.Krig, predict.se.Krig, predict.surface.se,
predict.surface, plot.Krig,
surface.Krig 
}
\examples{

# a 2-d example 
# fitting a surface to ozone  
# measurements. Exponential covariance, range parameter is 20 (in miles) 

fit <- Krig(ozone$x, ozone$y, theta=20)  
 
summary( fit) # summary of fit 
set.panel( 2,2) 
plot(fit) # four diagnostic plots of fit  
set.panel()
surface( fit, type="C") # look at the surface 

# predict at data
predict( fit)

# predict on a grid ( grid chosen here by defaults)
 out<- predict.surface( fit)
 surface( out, type="C") # option "C" our favorite

# predict at arbitrary points (10,-10) and (20, 15)
 xnew<- rbind( c( 10, -10), c( 20, 15))
 predict( fit, xnew)

# standard errors of prediction based on covariance model.  
 predict.se( fit, xnew)

# surface of standard errors on a default grid
 predict.surface.se( fit)-> out.p # this takes some time!
 surface( out.p, type="C")
 points( fit$x)


# Using anohter stationary covariance. 
# smoothness is the shape parameter for the Matern. 

fit <- Krig(ozone$x, ozone$y, Covariance="Matern", theta=10, smoothness=1.0)  
summary( fit)

#
# Roll your own: creating very simple user defined Gaussian covariance 
#

test.cov <- function(x1,x2,theta,marginal=FALSE,C=NA){
   # return marginal variance
     if( marginal) { return(rep( 1, nrow( x1)))}

    # find cross covariance matrix     
      temp<- exp(-(rdist(x1,x2)/theta)**2)
      if( is.na(C[1])){
          return( temp)}
      else{
          return( temp\%*\%C)}
      } 
#
# use this and put in quadratic polynomial fixed function 


 fit.flame<- Krig(flame$x, flame$y, cov.function="test.cov", m=3, theta=.5)

#
# note how range parameter is passed to Krig.   
# BTW:  GCV indicates an interpolating model (nugget variance is zero) 
#

# take a look ...
 surface(fit.flame, type="I") 

# 
# Thin plate spline fit to ozone data using the radial 
# basis function as a generalized covariance function 
#
# p=2 is the power in the radial basis function (with a log term added for 
# even dimensions)
# If m is the degree of derivative in penalty then p=2m-d 
# where d is the dimension of x. p must be greater than 0. 
#  In the example below p = 2*2 - 2 = 2  
#

 out<- Krig( ozone$x, ozone$y,cov.function="Rad.cov", 
                       m=2,p=2,scale.type="range") 

# See also the Fields function Tps
# out  should be identical to  Tps( ozone$x, ozone$y)
# 

# A Knot example

 data(ozone2)
 y16<- ozone2$y[16,] 

# there are some missing values -- remove them 
 good<- !is.na( y16)
 y<- y16[good] 
 x<- ozone2$lon.lat[ good,]

#
# the knots can be arbitrary but just for fun find them with a space 
# filling design. Here we select  50 from the full set of 147 points
#
 xknots<- cover.design( x, 50, num.nn= 75)$design  # select 50 knot points

 out<- Krig( x, y, knots=xknots,  cov.function="Exp.cov", theta=300)  
 summary( out)
# note that that trA found by GCV is around 17 so 50>17  knots may be a 
# reasonable approximation to the full estimator. 
#

# the plot 
 surface( out, type="C")
 US( add=TRUE)
 points( x, col=2)
 points( xknots, cex=2, pch="O")



# A correlation model example

# fit krig surface using a mean and sd function to standardize 
# first get stats from 1987 summer Midwest O3 data set 
# Compare the function Tps to the call to Krig given above 
# fit tps surfaces to the mean and sd  points.  
# (a shortcut is being taken here just using the lon/lat coordinates) 

 data(ozone2)
 stats.o3<- stats( ozone2$y)
 mean.o3<- Tps( ozone2$lon.lat, c( stats.o3[2,]))
 sd.o3<- Tps(  ozone2$lon.lat, c( stats.o3[3,]))

#
# Now use these to fit particular day ( day 16) 
# and use great circle distance
#NOTE: there are some missing values for day 16. 

 fit<- Krig( ozone2$lon.lat, y16, 
            theta=350, mean.obj=mean.o3, sd.obj=sd.o3, 
            Covariance="Matern", Distance="rdist.earth",
            smoothness=1.0,
            na.rm=TRUE) #


# the finale
 surface( fit, type="I")
 US( add=TRUE)
 points( fit$x)
 title("Estimated ozone surface")

#
#
# explore some different values for the range and lambda using REML
 theta <- seq( 300,400,,10)
 PLL<- matrix( NA, 10,80)
# the loop 
 for( k in 1:10){

# call to Krig with different ranges
# also turn off warnings for GCV search 
# to avoid lots of messages. (not recommended in general!)

  PLL[k,]<- Krig( ozone2$lon.lat[good,], y16[good],
             cov.function="stationary.cov", 
             theta=theta[k], mean.obj=mean.o3, sd.obj=sd.o3, 
             Covariance="Matern",smoothness=.5, 
             Distance="rdist.earth", nstep.cv=80,
             give.warnings=FALSE)$gcv.grid[,7]
  
#
# gcv.grid is the grid search output from 
# the optimization for estimating different estimates for lambda including 
# REML
# default grid is equally spaced in eff.df scale ( and should the same across theta)
#  here 

 }

# see the 2 column of $gcv.grid to get the effective degress of freedom. 

 cat( "all done", fill=TRUE)
 contour( theta, 1:80, PLL)
 

}

\keyword{spatial}
% docclass is function
