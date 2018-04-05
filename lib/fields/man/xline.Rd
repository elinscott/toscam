\name{xline}
\alias{xline}
\title{
  Draw a vertical line 
}
\description{
Adds  vertical lines in the plot region. 
}
\usage{
xline(x, ...)
}
\arguments{
\item{x}{
Values on x axis specifying location of vertical lines. 
}
\item{\dots}{
Any ploting options for abline. 
}
}
\seealso{
 yline, abline 
}
\examples{

plot( 1:10)
xline( 6.5, col=2)
 
world( col=3) 
yline( seq( -80,80,10),col=4, lty=2)
xline( seq( -180,180,10),col=4,lty=2)
yline( 0, lwd=2, col=4)
}
\keyword{aplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
