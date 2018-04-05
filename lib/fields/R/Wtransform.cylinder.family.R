"WQS.periodic" <-
function (x) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(-1, 3, 3, -1)
    D.rough <- c(-1, 3, -3, 1)
    D.le <- matrix(c(D.smooth, D.rough), nrow = 2, ncol= 4, 
        byrow = TRUE)
    D.re <- D.le
    n <- dim(x)[1]
    m <- dim(x)[2]
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:2, ] <- D.le %*% x[c(n, 1:3), ]
    tmp[(n - 1):n, ] <- D.re %*% x[c((n - 2):n,1), ]
    stuff <- n - 4
    indx <- seq(2, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp[c(seq(1, n, 2), seq(2, n, 2)), ]
}
"WQS.periodic.basis" <-
function (N, cut.n = 8) 
{
    x <- diag(1, N)
    nrow <- nrow(x)
    ncol <- ncol(x)
    nn <- nrow
    temp <- x
    nn <- cut.n * 2
    while (nn <= nrow) {
        temp[1:nn, ] <- WQSi.periodic(temp[1:nn, ])
        nn <- nn * 2
    }
    return(temp)
}
"WQS.periodic.T" <-
function (x) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <-  c(3, -3, -1, -1)
    D.rough <-  c(-1, 1, 3, 3)
    D.le <- matrix(D.smooth, nrow = 1, ncol = 4, byrow = TRUE)
    D.re <- matrix(D.rough, nrow = 1, ncol = 4, byrow = TRUE)
    n <- dim(x)[1]
    m <- dim(x)[2]
    x[c(seq(1, n, 2), seq(2, n, 2)), ] <- x
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[n, ] <- D.le %*% x[c(n-1,n,1:2), ]
    tmp[1, ] <- D.re %*% x[c((n - 1),n,1,2), ]
    stuff <- n - 3
    indx <- seq(1, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp
}
"WQS2d.cylinder" <-
function (x, transpose = FALSE, byX=T) 
{
if( byX){
    if (!transpose) 
        t(WQS(t(WQS.periodic(x))))
    else {
        WQS.periodic.T(t(WQS.T(t(x))))
    }
}
else{
   if (!transpose)
        t(WQS.periodic(t(WQS(x))))
    else {
        WQS.T(t(WQS.periodic.T(t(x))))
    }
}

}
"WQS2di.cylinder" <-
function (x, transpose = FALSE, byX=T) 
{
   if( byX){
 if (!transpose) 
        t(WQSi(t(WQSi.periodic(x))))
    else {
        WQSi.periodic.T(t(WQSi.T(t(x))))
    }
}
else{
 if (!transpose)
        t(WQSi.periodic(t(WQSi(x))))
    else {
        WQSi.T(t(WQSi.periodic.T(t(x))))
    }
}

}
"WQSi.periodic" <-
function (x) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(3, -3, 1, 1)/16
    D.rough <- c(1, -1, 3, 3)/16
    D.le <- matrix(D.smooth, nrow = 1, ncol = 4, byrow = TRUE)
    D.re <- matrix(D.rough, nrow = 1, ncol = 4, byrow = TRUE)
    n <- dim(x)[1]
    m <- dim(x)[2]
    x[c(seq(1, n, 2), seq(2, n, 2)), ] <- x
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[n, ] <- D.le %*% x[c(n-1,n,1:2), ]
    tmp[1, ] <- D.re %*% x[c((n - 1),n,1,2), ]
    stuff <- n - 3
    indx <- seq(1, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp
}
"WQSi.periodic.T" <-
function (x) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <-  c(1, 3, 3, 1)/16
    D.rough <-  c(1, 3, -3, -1)/16
    D.le <- matrix(c(D.smooth, D.rough), nrow = 2, ncol= 4, 
        byrow = TRUE)
    D.re <- D.le
    n <- dim(x)[1]
    m <- dim(x)[2]
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:2, ] <- D.le %*% x[c(n, 1:3), ]
    tmp[(n - 1):n, ] <- D.re %*% x[c((n - 2):n,1), ]
    stuff <- n - 4
    indx <- seq(2, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp[c(seq(1, n, 2), seq(2, n, 2)), ]
}
"Wtransform.cylinder.image" <-
function (x, inv = FALSE, transpose = FALSE, byX=TRUE, cut.min = 8) 
{
#
    if (transpose){ 
        inv <- !inv}

    n <- dim(x)[1]
    m <- dim(x)[2]
    if (n > m) {
        flip <- TRUE
        temp <- t(x)
        n <- dim(temp)[1]
        m <- dim(temp)[2]
#
# if the calculations are done on the transpose of the image then
# switch the choice of which dimension has the periodic basis
#
        byX <- !byX
    }
    else {
        flip <- FALSE
        temp <- x
    }

    if (n > m) 
        stop(" number of columns of x must >= to number of\nrows")
    nn <- n
    mm <- m
# test for dimensions "close" to dyadic (see help file)

           if( dyadic.2check( mm,nn,cut.min)==FALSE) 
                 {stop("error in length of column or row dimensions")}
    if (!inv) {
        while (nn > cut.min) {
            if (!transpose) {
                temp[1:nn, 1:mm] <- WQS2d.cylinder(temp[1:nn,1:mm], 
byX=byX)
            }
            else {
                temp[1:nn, 1:mm] <- WQS2di.cylinder(temp[1:nn, 1:mm], 
                  transpose = TRUE, byX=byX)
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
                temp[1:nn, 1:mm] <- WQS2di.cylinder(temp[1:nn, 1:mm], 
byX=byX)
            }
            else {
                temp[1:nn, 1:mm] <- WQS2d.cylinder(temp[1:nn, 1:mm], 
transpose = TRUE, byX=byX)
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
