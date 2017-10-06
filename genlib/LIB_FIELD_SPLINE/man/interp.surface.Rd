\name{interp.surface}
\alias{interp.surface}
\alias{interp.surface.grid}
\title{
  Fast bilinear interpolator from a grid.
}
\description{
Uses bilinear weights to interpolate values on a rectangular
grid to arbitrary locations or to another grid.
}
\usage{
interp.surface(obj, loc)
interp.surface.grid(obj, grid.list)

}
\arguments{
\item{obj}{
A list with components x,y, and z in the same style as used by contour,
persp, image etc. x and y are the X and Y grid values and z is a matrix
with the corresponding values of the surface
}
\item{loc}{
A matrix of (irregular) locations to interpolate. First column of loc
isthe X coordinates and second is the Y's. 
}
\item{grid.list}{ A list with components x and y 
describing the grid to interpolate}

}
\value{
 An vector of interpolated values.  NA are returned for
regions of the obj\$z that are NA and also for locations outside of the
range of the parent grid. 
}
\details{
Here is a brief explanation of
the interpolation:  Suppose that the location, (locx, locy)  lies in
between the first two grid points in both x an y. That is  locx is between
x1 and x2 and
locy is between y1 and y2.  Let ex= (l1-x1)/(x2-x1) ey= (l2-y1)/(y2-y1).
The
interpolant is

( 1-ex)(1-ey)*z11 + (1- ex)(ey)*z12 + ( ex)(1-ey)*z21 + ( ex)(ey)*z22  

Where the z's are the corresponding elements of the Z matrix. 

Note that bilinear interpolation can produce some artifacts related to
the grid and not reproduce higher behavior in the surface. For, example
the extrema of the interpolated surface will always be at the parent
grid locations.  There is nothing special about about interpolating to
another grid, this function just includes a \code{for} loop over one
dimension and a call to the function for irregular locations. It was
included in fields for convenience. since the grid format is so common. 

See also the akima package for fast interpolation from irrgeular locations. 

}
\seealso{
image.smooth, as.surface, as.image, image.plot, krig.image,Tps
}
\examples{
#
# evaluate an image at a finer grid
# 

data( lennon)
# create the surface object
obj<- list( x= 1:20, y=1:20, z= lennon[ 201:220, 201:220])

# sample at 50 equally spaced points 
temp<- seq( 1,20,,50)
make.surface.grid( list( temp,temp))-> loc
interp.surface( obj, loc)-> look
# take a look
image.plot( as.surface( loc, look))

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
