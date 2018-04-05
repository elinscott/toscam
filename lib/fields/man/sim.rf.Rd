\name{sim.rf}
\alias{sim.rf}
\title{
  Simulates a random field  
}
\description{
Simulates a random Gaussian field on a regular grid. 
}
\usage{
sim.rf(obj)
}
\arguments{
\item{obj}{
A covariance object that includes information about the covariance function 
and the grid for evaluation. Usually this created by a setup call to 
Exp.image.cov. (See details below.) 
}
\item{\dots}{
Additional arguments passed to a particular method.}

}
\value{
A matrix with the random field values 
}
\details{
This function takes an object that includes some preliminary calculations 
and so is more efficient for simulating more than one field from the same 
covariance. However, the algorithm using a 2-d FFT may not always work if 
the correlation scale is large (See the FIELDS manual for more details.)
The simple fix is increase the size of the domain so that the correlation 
sale becomes smaller relative to the extent of th domain. 

For a stationary model the covariance object has the components: 

names( obj) 
 "m"    "n"    "grid" "N"    "M"    "wght"   

.
where m and n are the number of grid points in x and y grid is a list with 
the grid point values for x and y N and M is the size of the larger grid 
that is used for simulation ( usually M= 2*m and N=2*n) to minimize 
periodic effects. wght is a matrix from the FFT of the covariance 
function.  The easiest way to create this object is to use for example 
Exp.image.cov with setup=T ( see below). 

The classic reference for this algorithm is 
Wood, A.T.A. and Chan, G. (1994).
    Simulation of Stationary Gaussian Processes in [0,1]d . Journal of
Computational and Graphical Statistics, 3, 409-432. 

}
\seealso{
Exp.image.cov, matern.image.cov  
}
\examples{
#Simulate a Gaussian random field with an exponential covariance function,  
#range parameter = 2.0 and the domain is  [0,5]X [0,5] evaluating the 
#field at a 100X100 grid.  
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100)) 
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
# Now simulate another ... 
look2<- sim.rf( obj)
# take a look 
set.panel(2,1)
 image.plot( grid$x, grid$y, look) 
 title("simulated gaussian field")
 image.plot( grid$x, grid$y, look2) 
 title("another (independent) realization ...")
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
