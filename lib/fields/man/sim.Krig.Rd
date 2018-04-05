\name{sim.Krig}
\alias{sim.Krig.standard}
\alias{sim.Krig.grid}
\title{Conditonal simulation of a spatial process}
\description{
Generates exact (or approximate) random draws from the conditional 
distribution of a spatial process given specific observations. This is a 
useful way to characterize the uncertainty in the predicted process from 
data. This is known as conditional simulation in geostatistics or 
generating an ensemble prediction in the geosciences. sim.Krig.grid can 
generate a conditional sample for a large regular grid but is restricted 
to stationary correlation functions. 
 }
\usage{
sim.Krig.standard(object, xp, M = 1, verbose = FALSE, sigma2 = NA, rho = NA)

sim.Krig.grid(object, grid.list = NA, M = 1, nx = 40, ny = 40, xy=c(1,2), verbose = 
FALSE, 
sigma2 = NA, rho = NA, extrap = FALSE)
} 
%- maybe also 'usage' for other objects documented here. 
\arguments{
  \item{object}{ A Krig object  }
  \item{xp}{ Locations where to evaluate the conditional process. }
  \item{M}{Number of draws from conditional distribution. }
  \item{verbose}{If true prints out intermediate information. }
  \item{sigma2}{User specified value for nugget variance or measurement 
         error. See Details below.}
  \item{rho}{User specified value for sill, or multiplier of spatial 
           covariance function. See Detials below.} 
  \item{grid.list}{ Grid information for evaluating the conditional 
surface as a grid.list.}
  \item{nx}{ Number of grid points in x.  }
  \item{ny}{ Number of grid points in y.}

 \item{xy}{ A two element vector giving the positions for the "X" and "Y"
variables for the surface. The positions refer to the columns of the  location
matrix used to define the multidimensional surface from the Krig object. This argument is
provided in lieu of generating the grid list. If a 4 dimensional surface
is fit to data then \code{ xy= c(2,4)} will evaluate a surface using the
second and fourth variables with  variables 1 and 3 fixed at their median
values. NOTE: this argument is ignored if a grid.list argument is
passed.  }

\item{extrap}{ If FALSE conditional process is not evaluated outside 
the convex hull of observations. } 

}

\details{
These functions generate samples from a conditional multivariate 
distribution that describes the uncertainty in the estimated spatial 
process under Gaussian assumptions. An important approximation throughout 
these functions is that all covariance parameters are fixed at their 
estimated or prescribed values. 

Given a spatial process  Z(x)= P(x) + h(x) observed at 

Y.k = P(x.k) + h(x.k) + e.k

where P(x) is a low order, fixed polynomial and h(x) a Gaussian spatial 
process.
With Y= Y.1, ..., Y.N,
the goal is to sample the conditional distribution of the process. 
 
[Z(x) | Y ]

For fixed a covariance this is just a multivariate normal sampling problem.
\code{sim.Krig.standard} samples this conditional process at the points 
\code{xp} and is exact for  fixed covariance parameters.
\code{sim.Krig.grid} also assumes fixed covariance parameters and does 
approxiamte sampling on a grid. 

The outline of the algorithm is 

0) Find the spatial prediction at the unobserved locations based on the actual data. Call this
 Z.hat(x). 

1) Generate an unconditional spatial 
process and from this process simluate synthetic observations. 

2) Use the spatial prediction model ( using the true covariance) to estimate the 
spatial process at unobserved locations. 

3) Find the difference between the simulated  process and its prediction based on synthetic 
observations. Call this e(x). 

4)  Z.hat(x) + e(x) is a draw from [Z(x) | Y ].

\code{sim.Krig.standard} follows this algorithm exactly. 

\code{sim.Krig.grid} evaluates the
conditional surface on grid and simulates the values of h(x) off the grid using bilinear 
interpolation of the four nearest grid points. Because of this approximation it is important to 
choose the grid to be fine relative to the spacing of the observations. The advantage of this
approximation is that one can consider conditional simulation for large grids -- beyond the 
size possible with exact methods. Here the method for simulation is circulant embedding and
so is restricted to correlation stationary fields. 

}
\value{

For \code{sim.Krig.standard} a matrix with columns indexed by the locations in \code{xp} and
\code{M} rows.

For \code{sim.Krig.grid} a list with arguments \code{x} and \code{y} defining the grid
 locations in the usual manner and \code{z}  contains the values of the simulated conditional
 field(s).    \code{z} is a three dimesional array where the first two indices are "x" and "y" and 
the third index is between 1 and M and indexes the simulated fields.

}


\author{Doug Nychka}
\seealso{ sim.rf, Krig}
\examples{
data( ozone2)

set.seed( 399)

# fit to day 16 from Midwest ozone data set.
Krig( ozone2$lon.lat, ozone2$y[16,], Covariance="Matern", 
theta=1.0,smoothness=1.0, na.rm=TRUE)-> out

# NOTE theta =1.0 is not the best choice but 
# allows the sim.rf circulant embedding algorithm to 
# work without increasing the domain.

#six missing data locations
 xp<-  ozone2$lon.lat[ is.na(ozone2$y[16,]),]

# 50 draws from process at xp given the data 
# this is an exact calculation
 sim.Krig.standard( out,xp, M=50)-> sim.out

# Compare: stats(sim.out)[3,] to  Exact: predict.se( out, xp)

# simulations on a grid
# NOTE this is approximate due to the bilinear interpolation
# for simulating the unconditional random field. 

sim.Krig.grid(out,M=5)-> sim.out

# take a look at the ensemble members. 

predict.surface( out, grid= list( x=sim.out$x, y=sim.out$y))-> look

zr<- c( 40, 200)

set.panel( 3,2)
image.plot( look, zlim=zr)
title("mean surface")

for ( k in 1:5){
image( sim.out$x, sim.out$y, sim.out$z[,,k], col=tim.colors(), zlim =zr)
}

}
\keyword{spatial}
% at least one, from doc/KEYWORDS
