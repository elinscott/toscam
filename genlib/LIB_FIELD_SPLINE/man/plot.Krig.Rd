\name{plot.Krig}
\alias{plot.Krig}
\alias{plot.sreg}
\title{
  Diagnostic and summary plots of a Kriging or spline object 
}
\description{
Plots a series of four diagnostic plots that summarize the fit. 
}
\usage{
plot.Krig(x, digits=4, which= 1:4,...)
plot.sreg(x, digits = 4, which = 1:4, ...)


}
\arguments{
\item{x}{ A Krig or an sreg object}

\item{digits}{
Number of significant digits for the RMSE label. 
}
\item{which}{
A vector specifying by number which of the four plots to draw.  
1:4 plots all four. 
}
\item{\dots}{
Optional graphics arguments to pass to each plot. 
}
}
\details{
 This function creates four summary plots of the Krig or sreg object. The
default is to put these on separate pages. However if the screen is
already divided in some other fashion the plots will just be added
according to that scheme. This option is useful to compare to compare
several different model fits. 

The first is a scatterplot of predicted value against observed.

The second plot is "standardized" residuals against predicted value.  
Here we mean that the residuals are divided by the GCV estimate for sigma
and multiplied by the square root of any weights that have been specified.
In the case of a "correlation model" the residuals are also divided by the
marginal standard deviation from this model.

The third plot are the values of the GCV function against the effective
degrees of freedom. When there are replicate points several versions of
the GCV function may be plotted.  GCV function is with respect to the
standardized data if a correlation model is specified. A vertical line
indicates the minimium found.

The fourth plot is a histogram of the standardized residuals. 
For sreg if multiple lambdas are given plotted are boxplots of the  
residuals for each fit. 
}
\seealso{
Krig, summary.Krig, Tps, set.panel
}
\examples{

fit<-Krig(ozone$x, ozone$y, theta=200)  
# fitting a surface to ozone  
# measurements 

set.panel( 2,2)
plot(fit)

fit<-sreg(rat.diet$t,rat.diet$con)

# fit rat data
set.panel(2,2)
plot(fit)       

set.panel(1,1) # reset graphics window. 

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
