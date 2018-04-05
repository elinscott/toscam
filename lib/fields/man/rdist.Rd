\name{rdist}
\alias{rdist}
\title{
  Euclidean distance matrix 
}
\description{
Given two sets of locations computes the Euclidean distance  matrix
among all pairings. 
}
\usage{
rdist(x1, x2)
}
\arguments{
\item{x1}{
Matrix of first set of locations where each row gives the coordinates of a
particular
point.
}
\item{x2}{
Matrix of second set of locations where each row gives the coordinates of
a particular point. If this is missing x1 is used. 
}
}
\value{
The distance matrix if nrow(x1)=m and nrow(
x2)=n then the returned matrix will be mXn. 
}
\details{
Let D be the mXn distance matrix.  The elements are 
the Euclidean distances between the all locations x1[i,] and x2[j,].

D.ij = sqrt(  sum.k (( x1[i,k] - x2[j,k]) **2 ).


}
\seealso{
  Exp.cov, rdist.earth 
}
\examples{
out<- rdist( ozone$x)
# out is a 20X20 matrix.
out2<- rdist( ozone$x[1:5,], ozone$x[11:20,])
#out2 is a 5X10 matrix
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
