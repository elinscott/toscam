\name{plot.Wimage}
\alias{plot.Wimage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots 2-d wavelet coefficents by level and type}
\description{
Produces a panel of images using the \code{split.screen} tools that organize the 
wavelet coefficients from a 
multiresolution by location, resolution level and type. 
}
\usage{
plot.Wimage(x, cut.min, graphics.reset = TRUE, common.range = FALSE, 
    color.table = tim.colors(128), Nlevel = NULL, with.lines = FALSE, 
    omd.width = 0.2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix of coefficients }
  \item{cut.min}{defines minimum level }
  \item{graphics.reset}{ 
If TRUE will reset device to orignal settings including closing the  split.screen
mode.}
  \item{common.range}{ If TRUE image plots will be on common color scale.  }
  \item{color.table}{Color table to be used for image plots. Default is Tim Hoar's 
favorite. }
  \item{Nlevel}{Number of levels to plot. Default is to plot all. }
  \item{with.lines}{If TRUE will add white outlines of pixels. Default is FALSE to 
prevent from tiny pixels from turning white!}
  \item{omd.width}{ Fration of device surface width devoted to the vertical color 
legend strips. Default usually accomodates most axis labels.}
  \item{\dots}{Other Graphical parameters to be passed to the \code{par} function.}
}
\details{
As with most complicated graphical figures you basically get what it draws although 
it should be easy to modify this function for customization. 
The \code{split.screen} set of graphical functions are used to divide up the 
plotting real estate into subplots of different sizes.  
The user can 
experiment with different outer margin space (omd.width)
Use cex.axis argument in the 
call to change the size of the numerals in the color strip. 

By setting grahics.reset to FALSE the function returns a matrix giving the 
split screen ids to reference each of the individual plots. Use the \code{screen}
function to move to given plot and use a high level plotting function to overlay 
information. (See example below.)
Use \code{close.screen( all=TRUE)} to turn off split.screen mode for subsequent 
and normal plotting. 
   }

\author{Doug Nychka}

\examples{
 
old.par<- par(no.readonly=TRUE) # these functions may leave the device in 
# with some funny defaults

#
#multiresolution of John Lennon 
data(lennon) 
Wtransform.image( lennon, cut.min=16)->look 
plot.Wimage( look, cut.min=16)

# adding information 
plot.Wimage( look, cut.min=16, Nlevel=3, graphics.reset=FALSE)-> plot.layout
# plot.layout here is a 4X3 matrix with the screen numbers

# move to the smooth coefficients plot

 screen( plot.layout[1,1])
 plot( c(.5,16.5), c( .5,16.5), type="n",axes=FALSE) 
 points( 8,8, cex=2, pch="+")

# NOTE: just points( 8,8) will not work here. This has to do with split.screen not 
# reseting the plotting pars correctly. 

 box(col=6, lwd=2)# just for fun 

close.screen( all=TRUE)
par( old.par) # reset  to old settings 


}
\keyword{hplot}
