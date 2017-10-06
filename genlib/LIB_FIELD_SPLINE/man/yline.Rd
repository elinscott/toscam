\name{yline}
\alias{yline}
\title{
  Draw horizontal lines 
}
\description{
Adds  horizontal lines in the plot region. 
}
\usage{
yline(y, ...)
}
\arguments{
\item{y}{
Values on y axis specifying location of vertical lines. 
}
\item{\dots}{
Any ploting options for abline. 
}
}
\seealso{
 xline, abline 
}
\examples{
world( col=3)
yline( seq( -80,80,10),col=4, lty=2)
xline( seq( -180,180,10),col=4,lty=2)
yline( 0, lwd=2, col=4)
}
\keyword{aplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
