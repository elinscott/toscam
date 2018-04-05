\name{bplot.xy}
\alias{bplot.xy}
\title{
  Boxplots for conditional distribution 
}
\description{
Draws boxplots for y by binning on x. This gives a coarse, but quick,  
representation
of the conditional distrubtion of [Y|X] in terms of boxplots.  
}
\usage{
bplot.xy(x, y, N = 10, breaks = pretty(x, N, eps.correct = 1),
         style = "tukey", outlier = TRUE, plot = TRUE, xaxt =
         "s", ...)
}
\arguments{
\item{x}{
Vector to use for bin membership
}
\item{y}{
Vector to use for constructing boxplot statistics. 
}
\item{N}{
Number of bins on x. Default is 10.
}
\item{breaks}{
Break points defining bin boundaries. These can be unequally spaced. 
}
\item{style}{
Type of boxplot default is "tukey". The other choice is "quantile" where
the whiskers are drawn to the 5 and 95 percentiles instead being based on
the inner fences.
}
\item{xaxt}{
 
Plotting parameter for x-axis generation. Default is to produce
an x-axis.
}
\item{outlier}{
If true outliers (points beyond outer fences) will be added to the plots.
}
\item{plot}{
If false just returns a list with the statistics used for plotting the
box plots.
}
\item{\dots }{
Any other optional arguments passed to the bplot.obj function there 
are quite a few that are useful, see the help file for bplot for details. 
}
}
\seealso{
bplot, draw.bplot
}
\examples{
# condition on swim times to see how run times vary
bplot.xy( minitri$swim, minitri$run, N=5)

# bivariate normal corr= .6
set.seed( 123)
x<-rnorm( 2000)
y<- .6*x +  sqrt( 1- .6**2)*rnorm( 200)
#
#
bplot.xy( x,y, breaks=seq( -3, 3,,15) ,xlim =c(-4,4), ylim =c(-4,4))
points( x,y, pch=".", col=3)
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
