\name{Wtransform}
\alias{Wtransform.image}
\alias{Wtransform}
\alias{Wtransform.cylinder.image}
\title{
  Quadratic W wavelet transform for 1-d vectors or  rectangular or 
cylindrical images 
}
\description{
Finds the forward or inverse discrete wavelet transform using the
quadratic W basis. 
}
\usage{
Wtransform(x, inv = FALSE, transpose = FALSE, cut.min = 8)

Wtransform.image(x, inv=FALSE, transpose=FALSE, cut.min=8)

Wtransform.cylinder.image(x, inv=FALSE, transpose=FALSE, 
byX=TRUE,cut.min=8)

}
\arguments{
\item{x}{
Matrix to be transformed. For Wtransform the number of rows must be 
dyadic or m=p*2\^L where p is lees than or equal to cut.min. A 1-d 
transformation is effected on each column of x. 
For 
the 2-d transforms both the row and the column dimensions must 
satisfy this condition. A 2-d transformation is applied to the image.  

}
\item{inv}{
If true will compute the inverse transform default is false  
}
\item{transpose}{
If true will compute the transpose of transform default is false 
}
\item{cut.min}{
Minimum level of transformation.  cut.min=8 means that the coarsest level
will consist of 8 scale functions for the 1-d case and 64=8X8 scale 
functions centered on an 8X8 grid for the 2-d case.  
}
\item{byX}{
If TRUE enforces periodic boundary conditions on horizontal coordinate 
of the image and if FALSE
applies condition to vertical coordinate. (See details.)
}
}
\value{
A matrix the same size as x.   
}
\details{
The wavelet transform can be thought as matrix multiplication 
W \%*\% vec(x)  where vec(x) is the matrix x  stacked by columns.  
The inverse transform is  inv(W) \%*\% vec(x) and transpose is  
 t(W) \%*\% vec( x). One can interpret the columns of inv(W) as basis 
functions 
and they follow the usual pattern of translations and dilations of mother 
and father wavelets. (See example below.)


Another way of thinking of the transformation is by recursion:  apply a
smoother and a "rougher" to a vector and taking every other value. Now
apply the same operation to the smooth results, now half the length of the
previous vector. At each step one reduced the vector by a half and cut.min
specifies the size when to stop. The function \code{WQS} performs one step
of the recursion with smoother weights (-1, 3, 3, -1) and roughening
weights (-1, 3, -3, 1) away from the edges. \code{WQS} actually works for
a matrix where each column is transformed in this manner. Boundary
adjustments are made to preserve the shape of the basis functions. (See
example below.) By convention the returned matrix has the smoothed results
in the first half of the columns and the rough results in the second half.
The discrete wavelet transform performs this operation recursively on the
smoothed results until the smoothed vector is less than a set size.
Because each step reduces the size of the vector by 2 it only makes sense 
to
apply this algorithm to vectors whose length is dyadic or the product 
of a
small integer and a dyadic (e.g. 96=3*32). The precise tests are done by 
\code{dyadic.check} and \code{dyadic.2check} and if n is the dimension the 
constraint it that n=p*2\^L where L is less or equal to than cut.min. 
The 
transform will result in p scale (father) basis functions if cut.min is 
equal to p. 


At the end of the day this
recursive algorithm defines a linear transformation of the original image
to something that we call the wavelet decomposition. This is the full
matrix W mentioned above. It is also possible to express multiplication of
inv(W) and transposes by a similar recursive scheme and related sets of
filter weights. (Try \code{ WQSi( WQS(x))} as a test.) The reader is
referred to \code{WQSi} \code{WQS.T} and \code{WQSi.T} for the filtering
primitives. Finally it should be noted that the particular wavelets chosen
here are not orthogonal by have nice smooth properties and the father
wavelet resembles a Gaussian density while the mother looks like its
derivative.  (See example below for some plots of the implied basis
functions.)
Note that \code{Wtransform} is "vectorized" so that with little extra 
overhead one 
can do transforms for many separate 1-d series with one call. 
In particular, \code{Wtransform(diag(1,128),inv=TRUE)} will neatly 
generate the W matrix given above. 

For two dimensions for the basic step one applies WQS to
the columns of the matrix and then to the rows.  e.g.  
\code{t(WQS(t(WQS(x)))}. This primitive step is implemented in WQS2d. The
final algorithm calls WQS2d or its variants repeatedly on a matrix that
decreases by a factor of two in size along each dimension at each 
iteration.  


If the cylinder variant is specified the transform uses periodic
boundaries on one of the coordinates. This is useful for data on zonal
section of a sphere where a constant line of latitude should be periodic.  
For byX=T the wavelet basis functions wrap on the x-axis when an image
plot is made. (See example below.) This convention may cause some
confusion because R experts will know that the image plot rotates the the
matrix so that the (1,1) element is at the lower left corner. Thus
enforcing periodicity along the X-coordinate of the image pertains to the
columns of the matrix used to represent the image in R. (Compare
\code{matrix( 1:10,2,5)} to \code{image( matrix( 1:10,2,5))}.

}
\section{References}{
Nychka,D. , C. Wikle and J.A. Royle, (2002)  Multi-resolution models for
nonstationary spatial covariance functions Statistical Modelling 2
315-332.

See also the original report on W matrices:

Man Kam Kwong and P. T. Peter Tang, "W-Matrices, Nonorthogonal 
Multiresolution Analysis, and Finite Signals of Arbitrary Length," 
Preprint MCS-P449-0794, September 1994 

\url{ftp://info.mcs.anl.gov/pub/tech_reports/reports/P449.ps.Z}
}

\seealso{
Wtransform.image, WQS, WQSi   
}
\examples{
# W transform of a simple function
x<- seq( 0,1,256)
y<- x*(1-x)**3 + ifelse( x>.3, x,0)

Wtransform( y)-> out

# plot the implied wavelet basis functions
ID<- diag( 1, 256)
WQS.basis( 256)-> Wbasis
set.panel(2,2)
matplot( 1:256, Wbasis[,1:8], type="l", lty=c(1,2), col=2)
title("Father")
matplot( 1:256, Wbasis[,9:16], type="l", lty=c(1,2), col=2)
title("Mother")
matplot( 1:256, Wbasis[,17:32], type="l", lty=c(1,2), col=2)
title("Mother scaled by 2")
matplot( 1:256, Wbasis[,33:64], type="l", lty=c(1,2), col=2)
title("Mother scaled by 4")

set.panel( 1,1)

# test that the transform works


# Precise definition of what the transform is doing in terms of 
# explicit matrix multiplication  all of 
# these should be machine zero
# Note that the direct matrix multiplications will be substantially slower
# for large vectors.
# y<- rnorm( 256)
# y<- y /sqrt(mean( y**2))


#sqrt(mean( c( Wtransform(y, inv=TRUE) - Wbasis%*%y)**2))
#sqrt(mean( c(Wtransform(y, inv=TRUE, transpose=TRUE) - t(Wbasis)%*%y)**2))
#sqrt(mean( c(Wtransform(y) - solve(Wbasis)%*%y)**2))
#sqrt(mean( c(Wtransform(y, transpose=TRUE) - t(solve(Wbasis))%*%y)**2))

#
## 2-d examples
#


# Wtransform of John Lennon image
data(lennon)
look<- Wtransform.image( lennon)
# 
### take a look: 
 image.plot( look)

 
# a diagonal detail basis function 
# we find this by just multipling W by a unit vector!

temp<- matrix(0, nrow=32, ncol=32) 
temp[8,5]<- 1 
look<- Wtransform.image( temp , inv=TRUE, cut.min=4)
image( look)
title("diagonal detail W-wavelet")

#just for fun: redo this example for all indices up to 8!
#
#set.panel( 8,8)
#par( mar=c(0,0,0,0))
#for (  k in (1:8)){
#for ( j in (1:8)){
#temp<- matrix( 0 , nx,ny)
#temp[k,j] <- 1
#Wtransform.image( temp, inv=T, cut.min=cut.min)-> look
#image( look, axes=FALSE, xlab="", ylab="")
#box()
#}
#}


# examine a basis function to see periodic condition enforces along 
# X axis of image. 
temp<- matrix(0, nrow=32, ncol=32)
temp[8,5]<- 1
image( Wtransform.cylinder.image( temp , inv=TRUE, cut.min=4))
# now along Y-axis
image( Wtransform.cylinder.image( temp , inv=TRUE, cut.min=4, byX=FALSE))

# reset panel 
set.panel( 1,1)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
