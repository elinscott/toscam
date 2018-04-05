\name{Wimage.cov}
\alias{Wimage.cov}
\alias{Wtransform.image.cov}
\alias{Wimage.sim}
\title{Functions for W-transform based covariance models }
\description{
These functions are designed for 2-d Gaussian random fields on a large regular grid. 
They require a sparse or block diagonal model for the covariances among the wavelet 
coefficients. Wimage.cov multiplies a vector with the implied covariance matrix and 
Wtransform.image.cov has similar functionality to Wimage.cov but is the 
correct form for fitting surfaces using \code{kirg.image}.
 Wimage.sim generates a random field with the implied covariance. }
\usage{
Wimage.cov(Y = NULL, ind = NULL, H.obj, find.row = FALSE)
Wtransform.image.cov(ind1, ind2=ind1, Y, cov.obj)
Wimage.sim( H.obj)
}
\arguments{
  \item{Y}{ The vector for the multiplication, if ind is missing this should a 
 matrix with the same dimensions as the grid locations. }
  \item{ind}{ If \code{Y} is not full grid \code{ind} gives the index
location for the Y subset. See details below for indexing conventions.}
  \item{H.obj}{A list with components need to describe the grid size and
components of H. See details below for a description and the example.  }
  \item{find.row}{ If true the row of the covariance given by ind is
returned.}
  \item{ind2}{This is the same as\code{ind} in \code{Wimage.cov}. }
  \item{ind1}{Indices giving subset of the full grid for the result of the
multipication. (See details below). Note that the default is that ind2 is
set to ind1.  }
  \item{cov.obj}{This is the same as H.obj above.}
}
  
\details{
Note that Wimage.cov only returns the product of a vector times the covariance 
matrix. This is usually all that is really needed and finesses the 
memory problems 
of dealing with very large covaraince matrices.  But see the last example for a 
loop to extract a submatrix of the covariance. 

In notation, if Sigma is the full covariance matrix among all grid points.
and Y is a vector then the intended calculation in R syntax is: 
\code{ Sigma[ind1, ind2]\%*\% Y}. 
This is accomplished by padding Y with zeroes to be the full vector/grid, 
doing the full multiplication of Sigma and then throwing away all 
elements but those in ind1. Note that for Wimage.cov, ind1 is assumed to 
be the full vector/grid and ind is equal to ind2 (when find.row is FALSE). 


The way 2-d wavelet basis functions are indexed can be confusing. See the
help files for Wimage.info and Wtransform.image for more background. For
these functions, ind can either be a single vector and the index refers to
the grid points in column stacked (or ``vec'') format. 

If ind has a two columns these refer to the row/column of the grid. 
Currently only one sparse representation for H is implemented, a
block/diagonal strategy. The covariance has the form Wi*H *H * Wi.t. 
Where Wi is mnXmn the matrix of wavelet basis functions with each column
being a basis evaluated on the gird points and in stacked column format.
Wi.t is its transpose. H is represented as partitioned matrix where H0 is
a full matrix and the remaining rows only have a diagonal term.  If ind0
denotes the indices for the elements of H0 then in R notation: 

H0 = H[ ind0, ind0] 

To simplify the coding the diagonal elements are represented as an mXn matrix
with those locations corresponding to H0 set to 1. 
Accordingly, in R notation 

H1= diag(H), H1[ind0] <- 1

With this representation, multiplying H by a vector u becomes

u[ind0]<- H0\%*\% c( u[ind0]) 

u*H1

Note that u is assumed to be an mXn matrix based on the grid dimensions. 

Under the assumptions the  grid ( or image) is mXn the components for H.obj are
\code{m} , \code{n},  \code{cut.min}, specifying the coarsest level of wavelet 
basis functions, H0, ind, and H1. 

By altering the multiplication of H  within these and changing the H.obj 
or cov.obj list one can easily create other sparse representations for 
the covariance. 

  }

\value{
For Wimage.cov and Wimage.sim an mXn matrix with the same extent as the grid 
defined by H.obj.

For Wtransform.image.cov a vector with the same length as \code{ind1} 
and at these grid locations. 

}
\author{Doug Nychka }

\seealso{Wtransform.image, Wimage.info, Exp.image.cov, matern.image.cov}
\examples{

# This example  creates a block/diagonal for H that is close to a 
# Matern covariance
# The example while lengthy is a practical run through of the process 
# of creating a reasonable sparse representation for H. 

M<- 16 # a 16X16 grid of locations
MM<- M**2
x<- 1:M
y<- 1:M
grid<- list( x=x, y=y)
cut.min<- 4
# H0 is father and mothers at first level
# with cut.min=4

# get indices associated with H0
I<-rep( 1:cut.min, cut.min)
J<- sort( I)
Wimage.s2i(I, J, flavor=0, level=0, cut.min=cut.min, m=M, n=M, mat=FALSE)-> ind0
Wimage.s2i(I, J, flavor=1, level=1, cut.min=cut.min, m=M, n=M, mat=FALSE)-> ind1
Wimage.s2i(I, J, flavor=2, level=1, cut.min=cut.min, m=M, n=M, mat=FALSE)-> ind2
Wimage.s2i(I, J, flavor=3, level=1, cut.min=cut.min, m=M, n=M, mat=FALSE)-> ind3
IND<- c( ind0, ind1, ind2, ind3)
NH<- length( IND)   # H0 will be NH by NH in size. 


# the Matern covariance function range=3, smoothness=1.
cov.obj<- stationary.image.cov( grid=grid, setup=TRUE, theta=3, smoothness=1)

# find elements of D0= H0**2  using the formula
#  D0 = W Sigma W.t  where W is the  2-d W transform ( the inverse of W.i above).
#  and Sigma is the Matern covariance
#
# the following looping algorithms are a way to avoiding explicitly creating
# gigantic covariance matrices and wavelet basis matrices associated with large 
# grids. Of course in this example the grid is small enough (16X16) that the 
# matrices could be formed explicitly 
#

D0<- matrix( 0, NH,NH)
# fill in D0
 for ( k in 1:NH){
 tmp<- matrix( 0, M,M)
 tmp[IND[k]] <- 1
 hold<- Wtransform.image( tmp,  transpose=TRUE, cut.min=cut.min)
 hold<- stationary.image.cov(Y=hold, cov.obj=cov.obj)
 hold<- Wtransform.image(hold, cut.min=cut.min)
 D0[k,] <- hold[IND]
}

# sqrt of D0
temp<- eigen( D0, symmetric=TRUE)

# next line should be 
#  H0<-temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)

H0<- temp$vectors \%*\% diag(sqrt(temp$values))\%*\% t(temp$vectors)


# find H1
H1<- matrix(0,M,M)
for ( k in 1:MM){
 tmp<- matrix( 0, M,M)
 tmp[k] <- 1
 hold<- Wtransform.image( tmp,  transpose=TRUE, cut.min=cut.min)
 hold<- matern.image.cov(Y=hold, cov.obj=cov.obj)
 hold<- Wtransform.image(hold, cut.min=cut.min)
 H1[k] <- sqrt(hold[k]) }
# remember that H1 has the H0 diagonal values set to 1. 
H1[IND] <- 1

#OK good to go.  Create the H.obj list
H.obj<- list( m=M, n=M, ind0=IND, cut.min=cut.min, H0=H0, H1=H1)
#
#
#


# mutliply the covariance matrix by a random vector
tmp<- matrix( rnorm( 256), 16,16)
Wimage.cov( tmp, H.obj=H.obj)-> look

# generate a random field
Wimage.sim(H.obj)-> look
image.plot( look)


#A check that this really works!
#find  the  135 =(9-1)*16 +7 == (7,9) row of covariance matrix
# we know what this should look like. 

Wimage.cov( ind= 135, H.obj=H.obj, find.row=TRUE)-> look
image.plot( x,y,look) # the covariance of the grid point (7,9) with 
                      # all other grid points -- a bump centered at (7,9)

#multiply a vector by just a subset of Sigma

ind<- sample( 1:256,50) # random set of 50 indices
Y<- rnorm( 50) # random vector of length 50 
Wimage.cov(Y, ind= ind, H.obj=H.obj)[ind]-> look

# look is now of length 50 -- as expected. 
# In R notation, this finding Sigma[ind, ind] %*% Y

# OK suppose you really need Sigma[ind,ind] 
# e.g. in order for   solve( Sigma[ind,ind], u)  
# here is a loop to do it.

Sigma<- matrix( NA, 50,50)
for ( k in 1:50){
# for each member of subset find row and subset to just the columns identified by ind 
Sigma[k,] <- c( Wimage.cov( ind= ind[k], H.obj=H.obj, find.row=TRUE)[ind])
}





}
\keyword{spatial }% at least one, from doc/KEYWORDS
