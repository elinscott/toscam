\name{quilt.plot}
\alias{quilt.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Image plot for irregular spatial data. }
\description{
Given a vector of z values associated with 2-d locations this 
function produces an 
image-like plot where the locations are discretized to a grid  and the z 
values are coded as a color level from a color scale. }
\usage{
quilt.plot(x, y, z, nrow = 64, ncol = 64, grid = NULL, 
                     add.legend=TRUE,add=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of the x coordinates of the locations -or- a
     a 2 column matrix of the x-y coordinates. }
  \item{y}{ A vector of the y coordinates -or- if the locations are passed 
in x the z vector  }
  \item{z}{ Values of the variable to be plotted.}
  \item{nrow}{ Number of grid boxes in x.}
  \item{ncol}{ Number of grid boxes in y. }
  \item{grid}{ A grid in the form of a \code{grid list}.  }
   \item{add.legend} { If TRUE a legend color strip is added}
  \item{add}{ If FALSE add to existing plot.}
  \item{\dots}{ arguments to be passed to the image.plot function }
}
\details{
This function combines the discretization to an image by the function 
\code{as.image} and is then graphed by \code{image.plot}.
Locations that fall into the same grid box will have their z values 
averaged. 

A similar function exists in the lattice package and produces good
looking plots. The advantage of this fields version is that it uses the
standard R graphics functions and is written in R code. Also, the
aggregation to average values for z values in the same grid box allows
for different choices of grids. If two locations are very close,
separating them could result in very small boxes. 

As always, legend placement is never completely automatic. Place the 
legend independently for more control, perhaps using \code{image.plot}
in tandem with \code{split.screen} or enlarging the plot margin 
See \code{help(image.plot)} for examples of this function and these 
strategies. }
\author{D.Nychka}
\seealso{ as.image, image.plot, lattice, persp, drape.plot } 
\examples{

data( ozone2)
# plot 16 day of ozone data set

quilt.plot( ozone2$lon.lat, ozone2$y[16,])
US( add=TRUE, col="grey", lwd=2)

#
# and ... if you are fussy 
# do it again 
# quilt.plot( ozone2$lon.lat, ozone2$y[16,],add=TRUE)
# to draw over the state boundaries. 
#

### adding a legend strip "by hand"
par( mar=c( 5,5,5,10)) # save some room for the legend
quilt.plot( ozone2$lon.lat, ozone2$y[16,], add.legend=FALSE)
image.plot(ozone2$lon.lat, ozone2$y[16,],legend.only=TRUE)

  }
}
\keyword{hplot}% at least one, from doc/KEYWORDS
