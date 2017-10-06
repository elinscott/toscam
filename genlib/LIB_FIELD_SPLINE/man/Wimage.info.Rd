\name{Wimage.info}
\alias{Wimage.info}
\alias{Wimage.i2s}
\alias{Wimage.s2i}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds key indices related to a 2-d multiresolution}
\description{
Functions for finding the indices and other information for a 2-d basis 
in a multiresolution wavelet decomposition. 
}
\usage{
Wimage.info(m = 128, n = m, cut.min = 4)

Wimage.i2s(ind, m, n, cut.min)

Wimage.s2i(i, j, level, flavor, m, n, cut.min, mat = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Nmber of rows of image (x) }
  \item{n}{number of columns of image (y) }
  \item{cut.min}{ The minimum number of father basis functions along one axis. }
  \item{ind}{indices for the basis functions can either be a vector or 2 
column matrix. }
 \item{i}{ Vector of row locations for basis functions. Or a list
 with components i, j, level and flavor. }
 \item{j}{ Vector of column locations for basis functions.}
 \item{level}{ Resolution level of basis functions.}
 \item{flavor}{ Type of basis function ( 0=S, 1=H, 2=V, 3=Di). See 
 details below. }
 \item{mat}{If TRUE returned index will be a matrix where the basis 
functions are indexed by a row and 
column. 
}
}
\details{

 The wavelet coefficents are computed efficiently as a single large
matrix/image but this format make it difficult to identify the indices
ofspecific kinds of basis functions.  The wavelet basis functions as found
with a tensor product multiresolution.  ( e.g. \code{ Wtransform.image})
are orgainzied in levels of resolution and type.  At the coarsest level of
resolution are father wavelets (dentoed by S for "smooth") roughly
centered on a rectanglar grid of locations (the size of the grid
controlled by \code{cut.min}). At this resolution are three mother
wavelets that capture horizontal (H), vertical (V) and diagonal( Di)
structure at this scale. These basis functions are also centered on a
lattice of grid locations. Subsequent levels only use the H,V, Di mother
wavelet templates. For example given a 32X32 image with cut.min=4 There
are three levels of resolution and 1024 basis functions total.  The
coarsest level is 16 S basis functions on a regular 4X4 grid, 16 H, 16 V
and 16 Di basis function also on 4X4 grid. The next level has 64 H, V and
Di basis functions on an 8X8 grid and the final level has 256 H, V and Di
basis functions on a 16X16 grid, giving a grand total of 1024 basis
functions.

Without some additional calculations it is possible to organize the 
results in 
different resolutions. This function provides the necessary indexing information to 
do this. See the function \code{plot.Wimag} as an exmaple of how this is used. 

Indexing the image can happen in 3 ways. 1)as a structure where one 
specifies the location (i,j) , level, and flavor. 
2) as the row and column indices of the image. 
3) as a single index if the image is stacked as a long vector. 
The functions Wimage.s2i and Wimage.i2s convert indices between these 
forms. 

}
\value{
  Wimage.info returns a list where all indices pertain to location or
subsets of the mXn matrix of wavelet coefficients.

  \item{m}{Number of rows in image}
  \item{n}{Number of columns of image}
  \item{cut.min}{cut.min}
  \item{S}{ A vector with 4 elements.  S[1]:S[2] range of rows for father wavelets
and S[3]:S[4] range of columns for father wavelets.
\item{H}{A matrix where each row is a resolution level and cloumns give s ranges of 
row and column indices. e.g.   H[k,1:4] gives the ranges for the rows 
and columns for the horizontal basis functions at level k. }
\item{V}{ A matrix in teh same format as H that gives subsets for the vertical 
basis 
functions.}
\item{Di}{ A matrix in the same format as H that gives subsets for the 
diagonal basis functions.}
\item{L}{ A column matrix where rows index resolutaion and column give the 
grid size for each level of resolution.} 
\item{Lmax}{ Total number of resolution levels.}
\item{offset.m}{A 3 column matrix that has the row offsets that indicate the 
begining of a block of coeffiecients.
Rows are levels of resolution and columns are H,V,Di.}

\item{offset.n}{A 3 column matrix that has the column offsets that indicate the 
begining of a block of coeffiecients.
Rows are levels of resolution and columns are H,V,Di.}
       }
  
}
\author{Doug Nychka }
\seealso{ \code{plot.Wimage}, \code{Wtransform.image} }

\examples{
#Find a basis function.
# For a basis on a 64X64 image with cut.min=8 find a 
# horizontal basis function at second level of resolution in the (4,4)
# position. ( There are 16X16 horizontal basis functions at the 2 nd level
#
Wimage.s2i( i=3,j=2, level=2, flavor=1, m=64, n=64, cut.min=8)-> ind
tmp<- matrix( 0, 64,64)
tmp[ind] <- 1
Wtransform.image( tmp, cut.min=8, inv=TRUE)-> look
image.plot( look)

# A check of Wimage.i2s
Wimage.i2s( ind,m=64, n=64, cut.min=8)
# should get i=3,j=2, level=2, flavor=1

# complete check of functions
 Wimage.i2s( 1:512, cut.min=8, m=16, n=32)-> look
 Wimage.s2i( look, cut.min=8, m=16, n=32, mat=FALSE)-> look2
 sum( look2 - (1:512)) # sum should be zero 


# W transform of John Lennon image

     data(lennon)
     m<- nrow(lennon)
     n<- ncol(lennon)

# add an grey strip down middle columns

     lennon[ , 128:132]<- 120

     look<- Wtransform.image( lennon, cut.min=8)

# get info 
     info<- Wimage.info( n, m, cut.min=8)


# Zero out all V basis functions coefficients, 

tmp<- look

Vind<- info$V
for (lev in 1: info$Lmax){
tmp[ Vind[lev,1]: Vind[lev,2], Vind[lev,3]: Vind[lev,4]]<- 0 
}


look2<- Wtransform.image( tmp, cut.min=8, inv=TRUE)

# take a look at vertically  filtered image
set.panel( 2,1)
par( pty="s")
image( lennon, col=grey(seq(0,1,,256)))
image( look2, col=grey(seq(0,1,,256)))

# NOTE:
# What is vertical and horizontal is confusing here
# due to the convention of the R image plot of running rows 
# of the image as the X coordinate
# The grey column in lennon image is plotted as 
# a horizontal line  
# 

}
\keyword{spatial}
