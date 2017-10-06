\name{tim.colors}
\alias{tim.colors}
\title{ Tim's useful color table }
\description{
A pleasing rainbow style color table patterned after that used in 
Matlab by Tim. 
}
\usage{
tim.colors(n = 64)
}
\arguments{
  \item{n}{ Number of color levels. The setting \code{n}=64 is the orignal 
definition.}
}
\details{
Ask Tim. 
}
\value{
A vector of character strings giving the colors in a hexadecimal format. 

}
\seealso{ topo.colors, terrain.colors }
\examples{

tim.colors(10) 
# returns an array of 10 strings in hex format
#e.g. (red, green,  blue) values of   (16,255, 239)
# translates to "#10FFEF" .

 

image( outer( 1:20,1:20,"+"), col=tim.colors( 75))
}
\keyword{ aplot}
