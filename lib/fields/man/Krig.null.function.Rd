\name{Krig.null.function}
\alias{Krig.null.function}
\title{Default function to create fixed matrix part of spatial process model.}
\description{
Constructs a matrix of terms representing a low order polynomial 
and binds additional columns due to covariates ( the Z matrix)
}
\usage{
Krig.null.function(x, Z = NULL, drop.Z = FALSE, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Spatial locations }
  \item{Z}{ Other covariates to be associated with each location.}
  \item{drop.Z}{If TRUE only the low order polynomial part is created.  }
  \item{m}{ The polynomial order is (m-1). }
}
\value{
A matrix where the first columns are the polynomial terms and the
following columns are from Z. 

}
\details{
This function can be modified to produce a different fixed part of the 
spatial model. The arguments x, Z and drop.Z are required but other arguments
can  be passed as part of a list in null.args  in the call to Krig. 
}
\author{Doug Nychka }
}
\seealso{ Krig  }
\keyword{ spatial}% at least one, from doc/KEYWORDS
