\name{krig.image}
\alias{krig.image}
\title{
  Spatial process estimate for large irregular 2-d dats sets. 
}
\description{
Computes the spatial predictions for large numbers of irregularly spaced
observations using the standard Kriging equations. The main approximation
is that the locations are discretized to a regular grid, but the field
need not be observed at all grid boxes.  

 
In Bayesian terms this
function computes the posterior mean for the field given the observations
under the usual Gaussian assumptions for the fields and observations.  
The solution is found by the iterative solution of a large linear system
using the conjugate gradient algorithm (CGA). Part of the calculations
rely on discretizing the spatial locations to a regular grid to make use
of the FFT for fast multiplication of a covariance matrix with a vector.
}
\usage{
krig.image(x, Y, cov.function, m=NULL, n=NULL, lambda=0, start=NULL, 
tol=1e-05, kmax=25,  cov.obj=NULL, grid=NULL, 
weights=rep(1, length(Y)), verbose=FALSE, conv.verbose=FALSE, expand=1, ...)
}
\arguments{
\item{x}{
A 2 column matrix of observed locations  
}
\item{Y}{
Values of observed field. Missing values are omitted from computation.  
}
\item{cov.function}{
An S function that multiplies the covariance matrix by a vector. Two that 
are part of FIELDS are Exp.image.cov ( Exponential and Gaussian) and  
W.image.cov ( W transform covariance model) 
}
\item{lambda}{
The value of the smoothing parameter. Should be nonnegative. See the notes 
below for more information about this parameter 
}
\item{m}{
Number of grid points in the x axis. Default is to use the length of  
grid\$x.  
}
\item{n}{
Number of grid points in the y axis. Default is to use the length of 
grid\$y. 
}
\item{cov.obj}{
A covariance object that contains information to be used by the  
covariance function  specified above.   
If this is not specified this object will be created within krig.image. 
}
\item{grid}{
A list with components x and y that specify the grid points in the  
x and y directions. The default is to use the number of point specified by 
m and n and use the ranges from the observed locations.  
}
\item{start}{
Starting values for omega2 in the iterative algorithm. Default is zero.   
}
\item{tol}{
Convergence tolerance for CGA. 
}
\item{kmax}{
Maximum number of iterations for CGA 
}
\item{weights}{
This vector is proportional to the reciprocal variance of the 
measurement errors. The default is a vector of ones.  
}
\item{verbose}{
If true all kinds of stuff is printed out! Default is false of course. 
}
\item{conv.verbose}{
If true the convergence criterion is printed out at each iteration of the  
CGA. The values are scaled as the criterion divided by the tolerance. So 
the algorithm terminates when the values are less than one.   
}
\item{expand}{
The amount the grid should be expanded beyond the range of the observed 
data. For example expand 1.1 will give a range that is 10 \% larger on 
each 
end.  
}
\item{\dots}{
Any extra arguments are considered as information for the covariance 
function and are used to create the covariance object.  
}
}
\value{
An list object of class krig.image.  
An explanation of some components: 

\item{call }{
The calling sequence   
}
\item{cov.function}{
A copy of the covariance S function  
}
\item{na.ind}{
logical indicating missing values in Y 
}
\item{xraw}{
Passed spatial locations having removed missing values 
}
\item{y}{
Observations having omitted missing values 
}
\item{N}{
Length of y 
}
\item{weights}{
passed weights having omitted missing cases.  
}
\item{lambda  }{
}
\item{grid}{
list with components x an y indicating grid for discretization 
}
\item{cov.obj}{
List object to use with cov.function 
}
\item{m}{
Number of grid point in x axis 
}
\item{n}{
Number of grid point in y axis 
}
\item{index}{
A two column matrix indicating the indices of the closest grid point to 
each observed location.  
}
\item{x}{
Observed locations discretized to nearest grid point 
}
\item{yM}{
Observed values but with a weighted average replacing multiple values
associated with the same grid point.  }
\item{xM}{
Discretized locations associated with yM
}
\item{weightsM}{
Weight vector
associated with YM.  }
\item{uniquerows}{
Logical indicating which rows of x
are unique.  
}
\item{shat.rep}{
Pooled standard
deviation among observations that fall within the same gird points
}
\item{indexM }{
A two column
matrix indicating the indices of the closest grid point to each observed
location, yM.  }
\item{qr.T}{
QR
decomposition of the matrix of constant and linear terms at xM
}
\item{multAx}{
The S function
that is used for matrix multiplication in the CGA.  }
\item{omega2}{
Parameter vector that describes the
spatial process part of the conditional mean.  }
\item{converge}{
CGA convergence information
}
\item{beta}{
Constant, and the
two linear parameters for the fixed linear part of the model
}
\item{delta}{
Covariance matrix
times delta give the spatial predictions.  }
\item{rhohat, rho}{
Conditional on lambda the MLE for
the parameter multiplying the covariance function.  }
\item{sigma2, shat.MLE}{
Conditional on lambda the MLE
for the parameter dividing the weight function.  }
\item{surface}{
A list giving the predicted surface at
the grid points.  }
\item{fitted.values}{
Predicted values at true
locations
}
}
\details{
 From a functional point of
view krig.image and supporting functions are similar to the class Krig.
The main difference is that only 2-dimensional problems are considered and
the solution is calculated for a fixed value of lambda. (The Krig function
can estimate lambda.)  For large data sets a practical way to estimate
lambda is by out of sample cross-validation and the FIELDS manual gives a
detailed example of this for the precip data set. Also see the manual for
an explanation of the computational strategy (Conjugate Gradient) here.

 Efficiency for large datasets comes with restrictions on the range of
covariance functions and some other features. Currently FIELDS just has
two covarince models: exponential/Gaussian and wavelet based. However, it
is not difficult to modify these to other models. The default
discretization is to a 64X64 grid however even 256X256 is manageable and
quite likely to separate irregular locations in most cases. The user
should also keep in mind that the estimate is the result of an iterative
algorithm and so issues such as good starting values and whether the
algorithm converged are present.

  
The spatial model includes a linear spatial drift and MLE estimates of the 
nugget variance and sill are found based on the values of lambda.  If the 
weights are all equal to one and the covariance function is actually a 
correlation function, in the notation of this function, the "sill" is 
sigma2 + rho and the "nugget" is sigma2. Moreover sigma2 and rho are 
constrained so sigma2/rho =lambda. This is why lambda is the crucial 
parameter in this model.  

Although the field is only estimated to the resolution of the grid, 
prediction off of the grid is supported by bilinear interpolation using
the FIELDS function interp.surface. 
}
\section{References}{
Large spatial prediction problems and nonstationary fields (1998) 
Nychka, D., Wikle, C. and Royle, J.A.  

FIELDS manual 
  
}
\seealso{
plot.krig.image, predict.krig.image, Exp.image.cov, sim.krig.image  
}
\examples{
# 
# fit a monthly precipitation  field over the Rocky Mountains 
# grid is 64X64 
out<- krig.image( x= RMprecip$x, Y = RMprecip$y, m=64,n=64,cov.function= 
Exp.image.cov,  
lambda=.5, theta=1, kmax=100)  

#  
# range parameter for exponential here is .5 degree in lon and lat.  
#diagnostic plots.  
plot( out)  

# look at the surface 
image.plot( out$surface)  #or just      surface( out)

#
#simulate 4 realizations from the conditional distribution 
look<- sim.krig.image( out, nreps=4) 
# take a look:  plot( look)
 
# check out another values of lambda reusing some of the objects from the 
# first fit 

out2<- krig.image( RMprecip$x, RMprecip$y, cov.function= Exp.image.cov, 
lambda=4, 
start= out$omega2,cov.obj=out$cov.obj)
#
# some of the obsare lumped together into a singel grid box
#
# find residuals among grid box means and predictions
res<- predict( out2, out2$xM) - out2$yM
#compare with sizes of out2$residuals (raw y data) 

#starting values from first fit  in out$omega2 
# covariance and grid information are 
# bundled in the cov.obj 
## 

#
## fitting a thin plate spline. The default here is a linear null space 
## and  second derivative type penalty term.  
##  you will just have to try different values of lambda vary them on  
##  log scale to 
 
out<- krig.image( RMprecip$x, RMprecip$y, cov.function=Rad.image.cov, 
lambda=1, m=64, n=64, p=2, kmax=300) 

# take a look: image.plot( out$surface) 

# check out  different values reuse some of the things to make it quicker 
# note addition of kmax argument to increase teh number of iterations

out2<- krig.image( RMprecip$x, RMprecip$y,cov.function=Rad.image.cov,  
lambda=.5, start= out$omega2, cov.obj=out$cov.obj, kmax=400) 

# here is something rougher  
out3<- krig.image( RMprecip$x, RMprecip$y,cov.function=Rad.image.cov, 
lambda=1e-2, start= out2$omega2, cov.obj=out$cov.obj,kmax=400,
tol=1e-3) 

# here is something close to an interpolation  
out4<- krig.image( RMprecip$x, RMprecip$y,cov.function=Rad.image.cov, 
lambda=1e-7, start= out3$omega2, cov.obj=out$cov.obj,kmax=500, tol=1e-3)
 
#compare the the four surfaces: 
# but note the differences in scales ( fix zlim to make them the same) 
# 
# take a look
# set.panel( 2,2) 
# image.plot( out$surface) 
# points( out$x, pch=".") 

# image.plot( out2$surface) 
# image.plot( out3$surface) 
# image.plot( out4$surface) 
 
# some diagnostic plots) 
set.panel( 4,4) 
plot( out, graphics.reset=FALSE) 
plot( out2, graphics.reset=FALSE) 
plot( out3, graphics.reset=FALSE) 
plot( out4, graphics.reset=FALSE) 
set.panel(1,1)
}

\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
