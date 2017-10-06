\name{vgram.matrix}
\alias{vgram.matrix}
\title{
  Computes a variogram from an image 
}
\description{
Computes a variogram for an image taking into account different directions 
and returning summary information about the differences in each of these 
directions.  
}
\usage{
vgram.matrix(dat, R=5, nsum=1:8, collapse=TRUE, dx = 1,dy = 1 )
}
\arguments{
\item{dat}{
A matrix spacing of rows and columns are assumed to have the same 
distance. 
}
\item{R}{
Maximum radius for finding variogram differences assuming that the grid
points are spaced one unit a part. Default is go out to a 
radius of 5.  
}
\item{nsum}{
The components returned from applying the fields stats function. 
Default is to return all 8 statistics, n, mean, s.d. quantiles etc.  
}
\item{collapse}{
If true the default just gives variogram as a function of distance between 
points.  
}
\item{dx}{
The spacing of grid points on the X axis. This is used to calculate the 
correct distance between grid points. If dx is not equal to dy then the 
collapse argument must be FALSE. 
}
\item{dy}{ The spacing of grid points on the Y axis.
See additional notes for dx.}

}
\value{
If collapse=T 
  
A list with components:  d, a vector of distances for the differences, 
and vgram, the variogram values. This is the traditional variogram 
ignoring direction.  
If collapse=F 
A list with components: d, a vector of distances for the differences, 
ind, a two column matrix giving the x and y increment used to compute the 
differences, stats, the result of applying the stats function to each set 
of differences,   and vgram, the mean of the differences 
}

\details{
When collapse=F this function is designed to look for departures from 
isotropy by separating the variogram differences according to 
orientation. For small R this runs efficiently because the differences are 
found by sub-setting the image matrix.  

For example, suppose that a row of the ind matrix is 
(2,3). The variogram value associated with this row is the mean of  
the differences (X(i,j)- X( i+2,j+3))**2 for all i and j. 
(Here X(.,.) are the values for the spatial field.)  In this example  
d= sqrt(13) and there will be another entry with the same distance 
but 
corresponding to the direction (3,2). 
plot.vgram.matrix attempts to organize all the different directions into a 
coherent image plot.
}
\seealso{
\code{\link{vgram}} 
}
\examples{
# variogram for Lennon image.
data(lennon)
out<-vgram.matrix( lennon) 

plot( out$d, out$vgram, xlab="separation distance", ylab="variogram") 
# image plot of vgram values by direction.  

# look at different directions 
out<-vgram.matrix( lennon, R=10, collapse=FALSE) #this takes a bit of time 

set.panel(2,1)
plot( out$d, out$vgram) 
plot(out$d, out$vgram.robust) 

#image plot of variogram values for different directions. 
set.panel(1,1)
plot.vgram.matrix( out)
# John Lennon appears remarkably isotropic!

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
