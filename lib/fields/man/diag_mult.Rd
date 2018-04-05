\name{\%d*\%-methods}
\docType{methods}
\alias{\%d*\%-methods}
\alias{\%d*\%,matrix,matrix-method}
\alias{\%d*\%,matrix,numeric-method}
\alias{\%d*\%,numeric,matrix-method}
\alias{\%d*\%,numeric,numeric-method}
\title{ Multiplying diagonal matrices  in Package `fields' }
\description{
Recognizes a left vector as the elements of a diagonal matrix and 
does the right multiplication efficiently (this is not rocket science!).
This method is used in the 
internal functions of Krig to make the code more readable. It avoids having 
a branch in the source code to handle the diagonal or nondiagonal cases.
Note that this operator is not symmetric: a vector in the left argument
is interpreted as a diagonal matrix and a vector in the right argument is 
kept as a column vector.  
}
\section{Methods}{
\describe{
\item{x = "matrix", y = "matrix"}{
computes x\%*\%y }

\item{x = "matrix", y = "numeric"}{
computes x\%*\% diag(y)}

\item{x = "numeric", y = "matrix"}{ 
computes diag( x)\%*\% y}

\item{x = "numeric", y = "numeric"}{ 
computes  diag(x)\%*\% y,  the diagonal matrix x multiplied by a vector. }

}

}
\keyword{methods}
