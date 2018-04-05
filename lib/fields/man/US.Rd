\name{US}
\alias{US}
\title{
  Plot of the US with state boundaries 
}
\description{
Plots quickly, medium resolution outlines of the US with the states and 
bodies of water.  
}
\usage{
US( xlim=c(-124.7, -67.1), ylim = c(25.2, 49.4), add=FALSE,shift=FALSE, ...)
}
\arguments{

\item{ylim}{
range of latitudes 
}
\item{xlim}{
range of longitudes 
}
\item{add}{
If true will add the world map to current plot 
}
\item{shift}{ If TRUE shifts to be centered on the Dateline and longitude
runs from 0 to 360. If FALSE centers on Prime Meridian and longitude runs
from -180 to 180.}

\item{\dots}{
These are graphical arguments that are passed to the lines function that 
draws outline.  
}
}
\details{
This function uses the FIELDS dataset US.dat for the coordinates.   
}
\seealso{
world
}
\examples{
# Draw map in device color # 3 
US( col=3) 
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
