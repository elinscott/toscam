\name{ozone}
\alias{ozone}
\title{
  Data set of ozone measurements at 20 Chicago monitoring stations. 
}
\description{
The \code{ozone} data is a list of 
components, x and y.  x 
component is longitude and latitude position of each of the 20 
Chicago monitoring stations, y is the average 
daily ozone values over the time period 6/3/87-8/30/87.
}

\format{
This data set is a list containing the following components: 

\describe{
 \item{lon.lat}{ Longitude-latitude positions of monitoring stations.  }
 \item{x}{An approximate Cartesian set of coordinates for the locations 
where the units are in miles. The origin is in the center of the 
locations. }
 \item{y}{ Average daily ozone values over 1987 summer. }
 }
}
\source{
AIRS, the EPA air quality data base. 
}
\seealso{
Tps, Krig  
}
\examples{
fit<- Tps(ozone$x, ozone$y) 
# fitting a surface to ozone measurements. 
surface( fit, type="I")
}
\keyword{datasets}
% docclass is function
% Converted by Sd2Rd version 1.21.
