"Wtransform" <-
function (x, inv = FALSE, transpose = FALSE, cut.min = 8) 
{
# coerce to one column matrix if x is not already a matrix
if( !is.matrix(x)){ x<- matrix( x, ncol=1)}

# tranpose operation requires similar recursion to the inverse
    if (transpose) 
        inv <- !inv

    nn<-n <- dim(x)[1]
    
#
# check that n is the product of a dyadic and an integer less or equal 
#to than cut.min
#
    if( dyadic.check( n,cut.min)==FALSE) 
{stop("error in column dimension ")}
    if (!inv) {
        while (nn > cut.min) {
            if (!transpose) {
                x[1:nn, ] <- WQS(x[1:nn, ])
            }
            else {
                x[1:nn, ] <- WQSi.T(x[1:nn, ]) 
            }
            nn <- nn/2
        }
    }
    if (inv) {
        NN <- n
        while (NN > cut.min) {
            NN <- NN/2
        }
        nn <- NN * 2

        while (nn <= n) {
            if (!transpose) {
                x[1:nn,] <- WQSi(x[1:nn, ])
            }
            else {
                x[1:nn,] <- WQS.T(x[1:nn, ])
            }
            nn <- nn * 2
        }
    }
return(x)
}
