\name{rdist.earth}
\alias{rdist.earth}
\title{
  Great circle distance matrix 
}
\description{
Given two sets of longitude/latitude locations computes the Great circle
(geographic)  distance matrix
among all pairings. 
}
\usage{
rdist.earth(x1, x2, miles = TRUE, R = NULL)
}
\arguments{
\item{x1}{
Matrix of first set of lon/lat coordinates first column is the
longitudes
and second is the latitudes. 
}
\item{x2}{
Matrix of second  set of lon/lat coordinates first column is the
longitudes
and second is the latitudes. If missing x1 is used. 
}
\item{miles}{
If true distances are in statute miles if false distances in kilometers. 
}
\item{R}{
Radius to use for sphere to find spherical distances. If NULL the radius
is either in miles or kilometers depending on the values of the miles
argument. If R=1 then distances are of course in radians. 
}
}
\value{
The great circle distance matrix if nrow(x1)=m and nrow(
x2)=n then the returned matrix will be mXn. 
}
\details{
Surprisingly this all done efficiently in R by dot products of the
direction cosines. Thanks to Qing Yang for pointing this out a long time 
ago. }
\seealso{
  rdist, stationary.cov 
}
\examples{
out<- rdist.earth ( ozone$lon.lat)
#out is a 20X20 distance matrix
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
