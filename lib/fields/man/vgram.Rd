\name{vgram}
\alias{vgram}
\title{
  Finds a traditional or robust variogram for spatial data. 
}
\description{
Computes pairwise squared differences as a function of distance. 
Returns either raw values or statistics from binning. 
}
\usage{
vgram(loc, y, id=NULL, d=NULL, lon.lat=FALSE, dmax=NULL, N=NULL, breaks=NULL)
}
\arguments{
\item{loc}{
Matrix where each row is the coordinates of an observed point
of the field
}
\item{y}{
Value of the field at locations
}
\item{id}{
A 2 column matrix that specifies which variogram differnces to find. 
If omitted all possible pairing are found. 
This can used if the data has an additional covariate that determines
proximity, for example a time window. 
}
\item{d}{
Distances among pairs indexed by id. If not included distances from from
directly from loc. 
}
\item{lon.lat }{
If true, locations are assumed to be longitudes and latitudes
and 
distances found are great circle distances ( in miles see rdist.earth).
Default is false. 
}
\item{dmax}{
Maximum distance to compute variogram.
}
\item{N}{
Number of bins to use.
}
\item{breaks}{
Bin boundaries for binning variogram values. Need not be equally spaced
but must be ordered. 
}
}
\value{

A list with these components. 

\item{vgram}
variogram values 

\item{d} pairwise distances 


\item{call}
calling string 

\item{stats}
Matrix of statistics for values in each bin. 
Rows are the summaries returned by the stats function or describe.  
If not either breaks or N arguments are not supplied then this
component is not computed. 

\item{centers}
Bin centers

}

\section{References}{
See any standard reference on spatial statistics. For example 
Cressie, Spatial Statistics  
}
\seealso{
vgram.matrix
bplot.xy, vgram.matrix 
}
\examples{
#
# compute variogram for the midwest ozone field  day 16
# (BTW this looks a bit strange!)
#
data( ozone2)
good<- !is.na(ozone2$y[16,])
x<- ozone2$lon.lat[good,] 
y<- ozone2$y[16,good]

look<-vgram( x,y, N=15, lon.lat=TRUE) # locations are in lon/lat so use right
#distance
# take a look:
#plot( look$d, look$vgram)
#lines(look$centers, look$stats["mean",], col=4)

brk<- seq( 0, 250,,25)
 
## or some boxplot bin summaries

bplot.xy( look$d, sqrt(look$vgram), breaks=brk,ylab="sqrt(VG)")
lines(look$centers, look$stats["mean",], col=4)

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
