\name{fields}
\alias{fields}
\title{
fields - tools for spatial data 
}
\description{
Fields is a collection of programs for curve and function
 fitting with an emphasis on spatial data and spatial statistics. The
 major methods implemented include cubic and thin plate splines, universal
 Kriging and Kriging for large data sets. One main feature is any
 covariance function implemented in R can be used for spatial prediction. 

fields stives to have readable and tutorial code. Take a look at the
source code for \code{Krig} and \code{Krig.engine.default} to see how
things work "under the hood". 

Some major methods include: 
\itemize{ \item \code{Tps} Thin Plate spline
regression (including GCV) 
\item \code{Krig} Spatial process estimation
(Kriging) including support for conditional simulation. 
}

The Krig function allow you to supply a covariance function that is
written in native S/R code. 
See (\code{stationary.cov}) that includes several families of 
covariances and distance metrics including the Matern and great circle
distance.  

Some other noteworthy functions are
\itemize{
\item \code{cover.design}  Space-filling designs where the distance 
function is expresed in R/S code

\item \code{as.image}, \code{image.plot}, \code{drape.plot}, \code{quilt.plot}
 
convenient functions for working with image data and rationally (well,
maybe reasonably) placing a color scale on an image plot. 

\item \code{sreg},  \code{qsreg} \code{splint}   Fast 1-D smoothing 
splines and 1-D 
quantile/robust and interpolating cubic splines
}

There are also generic functions that support 
these methods such as 

\code{plot} - diagnostic plots of fit \cr
\code{summary}- statistical summary of fit \cr
\code{print}- shorter version of summary \cr
\code{surface}- graphical display of fitted surface \cr
\code{predict}- evaluation fit at arbitrary points \cr
\code{predict.se}- prediction standard errors at arbitrary points. \cr
\code{sim.rf}- Simulate a random fields on a 2-d grid.

To get started, try some of the examples from help files for \code{Tps} or 
\code{Krig}. See also the manual/tutorial at 
\url{http://www.cgd.ucar.edu/stats/Software/fields}

Graphics tips:
\code{help( fields.hints)}
 gives some R code tricks for setting up common legends and axes. 
And has little to do with this package!

Testing:
See \code{help(fields.tests)} for testing fields. 

DISCLAIMER:

This is software for statistical research and not for commercial uses. The
authors do not guarantee the correctness of any function or program in
this package. Any changes to the software should not be made without the
authors permission.
}
\examples{
set.panel()
# some air quality data:
data(ozone2)
x<-ozone2$lon.lat
y<- ozone2$y[16,]

# pixel plot of spatial data
quilt.plot( x,y)
US( add=TRUE) # add US map

fit<- Tps(x,y)
# fits a GCV thin plate smoothing spline surface to ozone measurements.
# Hey, it does not get any easier than this!

summary(fit) #diagnostic summary of the fit 

set.panel(2,2)
plot(fit) # four diagnostic plots of  fit and residuals.

surface(fit) # contour/image plot of the fitted surface
US( add=TRUE, col="magenta", lwd=2) # US map overlaid
title("Daily max 8 hour ozone in PPB,  June 19th, 1987")
}
\keyword{datasets}
