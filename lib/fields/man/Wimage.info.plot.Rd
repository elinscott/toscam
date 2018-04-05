\name{Wimage.info.plot}
\alias{Wimage.info.plot}

\title{Plot to check 2-d multiresolution indexing}
\description{
Plots in image format the differnt blocks of coefficents associated 
with a 2-d multiresolution. This function can be used to check the 
Wtransform indexing 
functions and is also an introduction to how the coefficients are 
organized. 
}
\usage{
Wimage.info.plot(m, n, cut.min)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Number of rows in image}
  \item{n}{Number of columns of image }
  \item{cut.min}{Smallest number of smooth basis functions. Coarsest resolution will 
have atleast cut.min X cut.min basis functions with support on a regular grid. }
}
\details{
This function was used to check the (compicated) indexing functions work. 
But it also might a useful graphics to desribe how the different levels of wavelet 
coefficients are pack into the image format. 
}
\author{Doug Nychka}
\seealso{ 
\code{\link{Wtransform.image}}, 
\code{\link{W.info}}, 
\code{\link{Wimage.info}}, 
}
\examples{
#
# 64X 128 image coarsest level has 8X16 smooth basis funcitons. 
#
# NOTE as a matrix the image plot is upside down! Rows are along X-axis and 
# columns on Y
# e.g. the (1,1) element is the lower left corner. 
#

Wimage.info.plot( 64, 128, cut.min=8)

}
\keyword{spatial}% at least one, from doc/KEYWORDS
