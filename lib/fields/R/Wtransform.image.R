"Wtransform.image" <-
function (x, inv = FALSE, transpose = FALSE, cut.min = 8) 
{
    if (transpose) 
        inv <- !inv
    n <- dim(x)[1]
    m <- dim(x)[2]
    if (n > m) {
        flip <- TRUE
        temp <- t(x)
        n <- dim(temp)[1]
        m <- dim(temp)[2]
    }
    else {
        flip <- FALSE
        temp <- x
    }
    if (n > m) 
        stop(" number of columns of x must >= to number of\nrows")
    nn <- n
    mm <- m
# test 
           if( dyadic.2check( mm,nn,cut.min)==FALSE) 
                 {stop("error in column or row dimensions")}


    if (!inv) {
        while (nn > cut.min) {
            if (!transpose) {
                temp[1:nn, 1:mm] <- WQS2d(temp[1:nn, 1:mm])
            }
            else {
                temp[1:nn, 1:mm] <- WQS2di(temp[1:nn, 1:mm], 
                  transpose = TRUE)
            }
            nn <- nn/2
            mm <- mm/2
        }
    }
    if (inv) {
        NN <- n
        MM <- m
        while (NN > cut.min) {
            NN <- NN/2
            MM <- MM/2
        }
        nn <- NN * 2
        mm <- MM * 2
        while (nn <= n) {
            if (!transpose) {
                temp[1:nn, 1:mm] <- WQS2di(temp[1:nn, 1:mm])
            }
            else {
                temp[1:nn, 1:mm] <- WQS2d(temp[1:nn, 1:mm], transpose = TRUE)
            }
            nn <- nn * 2
            mm <- mm * 2
        }
    }
    if (flip) {
        return(t(temp))
    }
    else {
        return(temp)
    }
}
