\name{poisson.cov}
\alias{poisson.cov}
\title{
  Poisson spherical covariance function 
}
\description{
Given two sets of locations in lon/lat computes the cross covariance 
matrix for
the Poisson covariance among all pairings. 
}
\usage{
poisson.cov(x1, x2, eta = .2)
}
\arguments{
\item{x1}{
Matrix of first set of locations where each row gives the coordinates of a
particular
point. First column is longitudes and the second column is latitudes.
}
\item{x2}{
Matrix of second set of locations where each row gives the coordinates of
a particular point. If this is missing x1 is used. 
}
\item{eta}{
Range (or scale) parameter. Should be in the interval [0,1]
}
}
\value{
If nrow(x1)=m and nrow(x2)=n then the returned cross covariance matrix, 
will be mXn. 
}
\details{
This covariance is one of the few closed form covariances for the sphere
and also know as the Poisson kernel.  If x1 and x2 are matrices where
nrow(x1)=m and nrow(x2)=n then this function should return a mXn matrix
where the (i,j) element is the covariance between the locations x1[i,] and
x2[j,]. The covariance is found as
 
(1-eta**2)/(1 - 2 * eta * D.ij + eta**2)**(1.5)

where  D.ij is the great circle 
distance between  x1[i,] and x2[j,].
}
\seealso{
Krig, rdist.earth
}
\examples{
# plot of covariance

x<- make.surface.grid( list( x=seq( -180,180,,40), y= seq( -85,85,,40)))
x0<- matrix( c(0,0), ncol=2)
look<- poisson.cov( x,x0, eta=.5)
image.plot(as.surface(x,look))
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
