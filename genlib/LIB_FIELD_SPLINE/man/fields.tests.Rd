\name{fields tests}
\alias{fields tests}
\alias{test.for.zero}

\title{
Testing fields functions
}
\description{

Some of the basic methods in feilds can be tested by directly
implementing the linear algebra using matrix expressions. These
comparisons are done in the the R source code test files in the tests
subdirectory of fields.  The function \code{test.for.zero} is useful for
comparing the tests in a meaninful and documented way. 

}

\usage{
test.for.zero( xtest,xtrue,  tol= 1.0e-8, relative=TRUE, tag=NULL)
}

\arguments{

\item{xtest}{Vector of target values}

\item{xtrue}{Vector of reference values}

\item{tol}{Tolerance to judge whether the test passes.}

\item{relative}{If true a relative error comparison is used. (See
 details below.)}

\item{tag}{ A text string to be printed out with the test results as a
   reference}

}

\details{
The scripts in the tests subdirectory are 

\describe{
\item{Krig.test.R:}{Tests basic parts of the Krig and Tps functions including replicated
and weighted observations. }

\item{Krig.se.test.R:}{Tests computations of standard errors for the
Kriging estimate.}

\item{Krig.se.grid.test.R}{Tests approximate standard errors for the
Krig function found by Monte Carlo conditional simulation.} 

\item{Krig.test.W.R}{Tests predictions and A matrix  when an off diagonal 
observation weight matrix is used.}

\item{Krig.se.W.R}{Tests standard errors when an off diagonal 
observation weight matrix is used.}

}

To run the tests just attach the fields library and source the testing file. 
In the the fields source code these are in a subdirectory "tests".



 \code{test.for.zero} is used to print out the result for each
individual comparison.
Failed tests are potentially very bad and are reported with a
string beginning 

"FAILED test value = ... "

If the object test.for.zero.flag exists ( it can have any value),
all the tests that pass print  text beginning,

 "    PASSED  test at tolerance ..."

This startegy means that if all tests succeed nothing and the
\code{object test.for.zero.flag} does not exist then nothing is printed
in the test scripts.  This is option simplifies the output scripts for
running through the tests -- no news is good news. 


FORM OF COMPARISON:
The actual test done is the sum of absolute differnces:

test value = \code{ 
sum( abs(c(xtest) - c( xtrue) ) ) /denom}

Where \code{demon} is either  \code{ mean( abs(c(xtrue)))} for relative error
or 1.0 otherwise. 

Note the use of "c" here to stack any structure in xtest and xtrue into
a vector. 

}
 
\keyword{misc}
