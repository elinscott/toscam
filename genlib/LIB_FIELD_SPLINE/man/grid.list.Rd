\name{grid list}
\alias{grid list}
\alias{grid.list}
\alias{fields.x.to.grid}
\alias{parse.grid.list}
\alias{fields.convert.grid}
\alias{discretize.image}

\title{
  Grid list format for describing equally spaced grids and some related
functions. 
}
\description{
    
The object grid.list refers to a list that contains information for  
evaluating a function on a 2-dimensional grid of points. If a function  
has more than two independent variables then one also needs to specify  
the constant levels for the variables that are not being varied. 
This format is used in several places in  fields  for 
functions that evaluate function estimates and plot surfaces.}

\usage{
parse.grid.list( grid.list, order.variables="xy")
fields.x.to.grid(x,nx=80, ny=80, xy=c(1,2))
fields.convert.grid( midpoint.grid )
discretize.image( x, m = 64, n = 64,  grid = NULL,
        expand = c(1, 1), boundary.grid=FALSE )

}
\arguments{
\item{grid.list}{ No surprises here -- a grid list! These can be unequally 
spaced.}
\item{order.variables}{ If "xy" the x variable will be subsequently plotted
as the horizontal variable. If "yx" the x variable will be on the vertical 
axis.}
\item{x} { A matrix of independent variables such as the locations of 
observations given to Krig.}
\item{nx}{Number of grid points for x variable.}
\item{ny}{Number of grid points for y variable.}
\item{m}{Number of grid points for x variable.}
\item{n}{Number of grid points for y variable.}
\item{xy}{The column positions that locate the x and y variables for the 
grid.} 
\item{grid}{ A grid list!}
\item{expand}{ 
A scalar or two column vector that will expand the
grid beyond the range of the observations.}

\item{midpoint.grid}{ Grid midpoints to convert to grid boundaries.}

\item{boundary.grid}{ If TRUE interpret grid points as boundaries  of 
grid boxes. If FALSE interpret as the midpoints of the boxes. }


}
\details{
The form of a grid.list is  
list( var.name1= what1  , var.name2=what2 , ... var.nameN=what3)  
Here var.names are the names of the independent variables. 
The what options describe what should be done with this variable when  
generating the grid. These should either an increasing sequence of points
or a single vaules. Obviously there should be only be two variables with 
sequences to define a grid for a surface.

Most of time the gridding  sequences  are equally  
spaced and are easily generated using the seq function. Also throughout fields
the grid points are typically the midpoints of the grid rather the grid box
boundaries. However, these functions can handle unequally spaced grids and the
logical boundary.grid can indicate a grid being the box boundaries. 

The variables in the list components are assumed  
to be in the same order as they appear in the data matrix. 


Some fields internal functions that support interpreting grid list format are:

\code{fields.x.to.grid}: 
 Takes an "x" matrix of locations or independent variables and creates a
reasonable grid list. This is used to evaluate predicted surfaces when a
grid list is not explicited given to predict.surface. The variables
(i.e. columns of x) that are not part of the grid are set to the median
values.  The x grid values are \code{nx} equally spaced points in the
range \code{x[, xy[1]]}. The y grid values are \code{ny} equally spaced
points in the range \code{x[, xy[2]]}.


\code{parse.grid.list}:
 Takes a grid list and returns the information in a more expanded list
form that is easy to use. This is used, for example, by predict.surface
to figure out what to do!


\code{fields.convert.grid}:
 Takes a vector of n values assumed to be midpoints of a grid and 
returns the n+1 boundaries. See how this is used in discretize.image 
with the cut function. This function will handle unequally spaced
grid values. 

 \code{discretize.image}: Takes a vector of locations and a 2-d grid and
figures out to which boxes they belong. The output matrix ind has the
grid locations. If boundary.grid is FALSE then the grid list (grid) is
assumed to be grid midpoints. The grid boundaries are taken to be the
point half way between these midpoints. The first and last boundaries
points are determined by extrapolating so that the first and last box
has the midpoint in its center. (See the code in fields.convert.grid for
details.) If grid is NULL then midpoints are found from m and n and the
range of the x matrix. 

}

\seealso{
as.surface, predict.surface, plot.surface, surface, 
expand.grid, as.image 
}

\examples{

#Given below are some examples of grid.list objects and the results  
#when they are used with make.surface.grid. Note that  
#make.surface.grid returns a matrix that retains the grid.list  
#information as an attribute. 

grid.l<- list( 1:3, 2:5) 
make.surface.grid(grid.l)

  
grid.l <- list( 1:3, 10, 1:3) 
make.surface.grid(grid.l) 

#The next  example shows how the grid.list can be used to  
#control surface plotting and evaluation of an estimated function. 
# first create a test function  

set.seed( 124)

X<- 2*cbind( runif(30), runif(30), runif(30)) -1
  
dimnames( X)<- list(NULL, c("X1","X2","X3")) 
y<- X[,1]**2 + X[,2]**2 + exp(X[,3])   

# fit an  interpolating thin plate spline  
out<- Tps( X,y) 

grid.l<- list( X1= seq( 0,1,,20), X2=.5, X3=seq(0,1,,25)) 
surface( out, grid.list=grid.l) 
#  surface plot based on a 20X25 grid in X1 an X3  
#                       over the square [0,2] and [0,2]   
#                       holding X2 equal to 1.0. 
#

discretize.image( RMprecip$x, m=15, n=15)-> look
Z<- matrix( 0, 15,15)

Z[look$ind]<- 1

image( look$grid$x, look$grid$y, Z) # indicator image of discretized locations. 
points( RMprecip$x,col="magenta") # actual locations
# (there may be more than one location in the grid boxes)

}
\keyword{misc}
% docclass is function
% Converted by Sd2Rd version 1.21.
