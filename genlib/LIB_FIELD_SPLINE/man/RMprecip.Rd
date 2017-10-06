\name{RMprecip}
\alias{RMprecip}
\title{
Monthly total precipitation (mm) for August 1963 in the Rocky Mountain 
Region 
}
\description{
This is a useful spatial data set of moderate size consisting of 865 
locations. See www.image.ucar.edu/stats/GSP/DATA for the source of these data. 

}

\format{
This data set is a list containing the following components: 

\describe{
\item{x}{
Longitude-latitude position of monitoring stations. Rows names are station id codes. 
}
\item{elev}{
Station elevation in meters. 
}
\item{y}{
Monthly total precipitation in millimeters.
for August 1963
}
}
}

\examples{
# this data set was created  the 
# historical data  taken from 
# Observed monthly precipitation, min and max temperatures for the coterminous US 
# 1895-1997
# NCAR_pinfill 
# see the Geophysical Statistics Project datasets page for the supporting functions 
# and details. 

# plot 
quilt.plot(RMprecip$x, RMprecip$y)
US( add=TRUE, col=2, lty=2)

#

} 
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
