\name{as.image}
\alias{as.image}
\title{
  Creates image from irregular x,y,z  
}
\description{
Discretizes a set of  2-d locations to a grid and produces a image object  
with the z values in the right cells. For cells with more than one Z
value the average is used.   
}
\usage{
as.image(Z, ind=NULL, grid=NULL, x=NULL, nrow=64, ncol=64,weights=NULL,
 na.rm=FALSE, nx=NULL,ny=NULL, boundary.grid=FALSE)
}
\arguments{
\item{Z}{
Values of image 
}
\item{ind}{
A matrix giving the row and column subscripts for each image 
value in Z. (Not needed if x is specified.) 
}
\item{grid}{
A list with components x and y of equally spaced values describing the 
centers of the grid points. The default is to use nrow and ncol and the 
ranges of the data locations (x)  to construct a grid.  
}
\item{x}{
Locations of image values. Not needed if  ind  is specified.  
}
\item{nrow}{
Number of rows in image matrix ( x-axis direction) 
}
\item{ncol}{
Number of columns in image matrix ( y-axis direction) 
}
\item{weights}{
If two or more values fall into the same  
pixel a weighted average is used to represent the pixel value. Default is 
equal weights. 
}
\item{na.rm}{ 
If true NA's are removed from the Z vector.}
\item{nx}{ Same as nrow}
\item{ny}{ Same as ncol}
\item{boundary.grid} { If FALSE grid points are assumed to be the 
midpoints. If TRUE they are the grid box boundaries.}
}
\value{
An list in image format with a few more components. Components x and y are 
the grid values , z is a 
nrow X ncol matrix 
with the Z values. NA's are placed at cell locations where Z data has 
not been supplied.  
Component ind is a 2 column matrix with subscripts for the locations of 
the values in the image matrix.  
Component N is an image matrix with the number of original values that 
fall within each grid box. 
Component weights is an image matrix with the  sum of the
individual weights for each cell.  
}
\details{
The discretization is straightforward once the grid is determined.  
If two or more Z values have locations in the same cell the average value 
is 
taken as the value.  See the source code to modify this to get more 
information about coincident locations. (See the call to fast.1way) 
}
\seealso{
image.smooth, image.plot, Krig.discretize, Krig.replicates  
}
\examples{
# convert precip data to 50X50 image  
look<- as.image( RMprecip$y, x= RMprecip$x, nrow=50, ncol=50)
image.plot( look) 

# number of obs in each cell.
image.plot( look$x ,look$y, look$N, col=terrain.colors(50)) 
# hot spot is around Denver
}
\keyword{manip}
% docclass is function
% Converted by Sd2Rd version 1.21.
