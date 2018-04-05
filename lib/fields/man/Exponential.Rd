\name{Exponential}
\alias{Exponential}
\alias{Matern}
\title{Covariance functions}
\description{
Functional form of covariance function assuming the argument is a
distance between locations. 
}
\usage{
Exponential(d, range = 1, alpha = 1/range, phi = 1)
Matern (d , scale = 1, range = 1,alpha=1/range,
     smoothness = 0.5, nu= smoothness, phi=scale) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{d}{ Vector of distances }

  \item{range}{ Range parameter default is one. Note 
       that the scale can also be specified  through the "theta" 
           scaling argument  used in fields covariance functions) }

  \item{alpha}{1/range }

  \item{scale}{Same as phi }

  \item{phi}{ Marginal variance.}

  \item{smoothness}{ Smoothness parameter in Matern. Controls the number
of derivatives in the process. Default is 1/2 corresponding to an exponential 
covariance.}

  \item{nu}{ Same as smoothness}

}

\details{
Exponential: 

phi* exp( -d/range)

Matern:

   phi*con*(d\^nu) * besselK(d , nu )

 Matern covariance function transcribed from Stein's book page 31
nu==smoothness, alpha ==  1/range

 GeoR parameters map to kappa==smoothness and phi == range
check for negative distances

\code{con} is a constant that normalizes the expression to be 1.0  when phi=1.0
and d=0. 

}
\value{

A vector of covariances. 

}
\references{ Stein's book}
\author{Doug Nychka}
\seealso{stationary.cov, stationary.image.cov}
}
\keyword{spatial}% at least one, from doc/KEYWORDS
