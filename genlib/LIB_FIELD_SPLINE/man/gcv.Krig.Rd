\name{gcv.Krig}
\alias{gcv.Krig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds profile likelihood and GCV estimates of 
smoothing parameters}
\description{
This is a secondary function that will use the computed Krig 
object and find various estimates of the 
smoothing parameter lambda. These are several different flavors of cross-validation,
a moment matching strategy  and
the profile likelihood. 
This function can also be used independently
with different data sets  (the y's) if the covariates ( the x's) are the same and
thus reduce the computation. 
}
\usage{
gcv.Krig(out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, 
                  verbose = FALSE, tol = 1e-05, offset = 0, 
                  y = NULL, give.warnings = TRUE, give.warnings.REML = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{ A Krig object}
  \item{lambda.grid}{ Grid of lambdas for coarse search. The default is 
equally spaced on effective degree of freedom scale. }
  \item{cost}{ Cost used in GCV denominator }
  \item{nstep.cv}{ Number of gird points in coarse search. }
  \item{rmse}{ Target root mean squared error to mathc with the estimate of sigma**2 }
  \item{verbose}{ If true prints intermediate results.  }
  \item{tol}{ Tolerance in delcaring convergence of golden section search or bisection search. }
  \item{offset}{ Additional degrees of freedom to be added into the GCV denominator.
}
  \item{y}{A new data vector to be used in place of the one associated with the
Krig object (obj) }
  \item{give.warnings}{ If FALSE will suppress warnings about grid search being out of 
range for various estimates based on GCV.}
  \item{give.warnings.REML}{ If FALSE will suppress warnings about grid search being 
out of 
range when finding REML estimate of lambda.}

}
\details{

 This function finds several estimates of the smoothing parameter using
first a coarse grid search followed by a refinement using a minimization (
in the case of GCV or maximum likelihood) or bisection in the case of
mathcing the rmse. Details of the estimators can be found in the help file
for the Krig function. 

The Krig object passed to this function has some matrix decompostions that
facilitate rapid computation of the GCV and ML functions and do not depend
on the independent variable. This makes it possible to compute the Krig
object once and to reuse the decompostions for multiple data sets. (But
keep in mind if the x values change then the object must be recalculated.)
The example below show show this can be used for a simulation study on the
variability for estimating the smoothing parameter. 


}
\value{A list giving a summary of estimates and diagonostic details with the 
following components:
  \item{gcv.grid }{ A matrix describing results of the 
coarse search rows are values of lambda and the columns are 
lambda= value of smoothing parameter, 
trA=effective degrees of freedom, 
GCV=Usual GCV criterion, 
GCV.one=GCV criterion leave-one-out, 
GCV.model= GCV based on average response in the case of replicates, 
shat= Implied estimate of sigma ,
-Log Profile= negative log of profiel likelihood for the lambda. 
}
  \item{lambda.est}{Summary table of all estimates
Rows index different types of estimates: 
GCV, GCV.model, GCV.one, RMSE, pure error, -Log Profile
and the columns are the estimated values  for lambda, trA, GCV, shat.
}
}
\author{Doug Nychka}

\seealso{ 
\code{\link{Krig}},
\code{\link{Tps}},
\code{\link{predict.Krig}} }
\examples{

# 
Tps( ozone$x, ozone$y)-> obj # default is to find lambda by GCV
summary( obj)

gcv.Krig( obj)-> out
print( out$lambda.est) # results agree with Tps summary

# a simulation example
x<- seq( 0,1,,150)
f<-  x**2*( 1-x)
f<- f/sqrt( var( f))

set.seed(123) # let's all use the same seed
sigma<- .1
y<- f + rnorm( 150)*sigma

Tps( x,y)-> obj # create Krig object

hold<- matrix( NA, ncol=4, nrow=20)

for( k in 1:20){
# look at GCV estimates of lambda
# new data simulated
   y<- f + rnorm(150)*sigma 
# save GCV estimates
hold[k,]<-  gcv.Krig(obj,  y=y, give.warnings=FALSE)$lambda.est[1,]
}
plot( hold[,2], hold[,4], xlab="estimated eff. df", ylab="sigma hat")
yline( sigma, col=2)
# note some occaisional flaky behaviour with GCV ( eff df > 20!)

  }

}
\keyword{spatial}
