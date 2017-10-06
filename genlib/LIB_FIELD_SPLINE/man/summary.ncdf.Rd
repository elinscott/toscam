\name{summary.ncdf}
\alias{summary.ncdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizes a netCDF file handle}
\description{
Provides a summary of the variable names and sizes from the handle
returned from netCDF file. 
}
\usage{
summary.ncdf(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
The "handle" returned by the \code{read.ncdf} function from the
ncdf package.  
 }
\item{\dots}{ Other arguments to pass to this function. Currently, no 
other arguments are used. 
}
}
\details{
This function is out of place in fields but was included because often large
geophysical data sets are in netCDF format and the ncdf R package is also 
needed. To date the summary capability in the ncdf package is limited and 
this function is used to supplement it use. The function is also a a 
useful device to see how the ncdf object is structured.
}

\author{ D. Nychka }

\seealso{ ncdf}
\keyword{ IO }% at least one, from doc/KEYWORDS
