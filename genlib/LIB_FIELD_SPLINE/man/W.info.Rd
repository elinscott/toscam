\name{W.info}
\alias{W.info}
\alias{W.i2s}
\alias{W.s2i}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gives indexing imfomration for a wavelet decompostion }
\description{
Functions for finding various indices and sizes of different 
parts of a 1-d multiresolution and converting to a sequential index.
}
\usage{
W.info(m = 128, cut.min = 4)

W.i2s(ind, m, cut.min)

W.s2i(i, level, m, cut.min)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Length of series }
  \item{cut.min}{Number of scale basis functions (or father wavelets) at the 
coarsest resolution ) }
  \item{ind}{Vector of indices of the basis functions.}
  \item{i}{ Position of the basis function within its level.}
  \item{level}{ level of resolution (father wavelets have level==0). }
 }
\details{
W.info gives summary information about the multiresolution. 
W.i2s converts a sequential index for the basis functions into a list with 
locations and levels. W.s2i is the inverse transformation from a list 
structure to an integer index. }
\value{
Return list for W.info
  \item{m }{length of series }
  \item{cut.min}{cut.min (what else!)}
  \item{S}{first and last indices for father wavelets}
  \item{H}{a matrix where rows are levels of resolution and columns are the first 
and last indices for the mother basis functions.}
  \item{L}{number of basis funcion in each level of the mother wavelets.}
   \item{Lmax} {the number of levels of resolution}
   \item{offset}{a matrix where rows are levels of resolution and the columns is the 
offset index for the beginning of the indices for basis funtion at that resolution 
level.
}   

}
\author{Doug Nychka}

\seealso{\code{Wtransform}, 
\code{Wtransform.image}, \code{plot.Wimage}, \code{Wimage.info}}
\examples{

# series of length 64 where the coarsest level is 8 father basis functions. 
W.info(64, 8)

#index for 4th basis location at the 2nd level

W.s2i( 4, 2, m=64, cut.min=8)
# location and level for the 48th basis function.
W.i2s( 48, m=64, cut.min=8)
}

\keyword{spatial}% at least one, from doc/KEYWORDS
