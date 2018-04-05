"WQS.T" <-
function (x) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(3, -3, -1, -1)
    D.le <- matrix(c(2, 2), nrow = 1, ncol = 2, byrow = TRUE)
    D.re <- matrix(c(2, -2), nrow = 1, ncol = 2, byrow = TRUE)
    D.rough <- c(-1, 1, 3, 3)
    n <- dim(x)[1]
    m <- dim(x)[2]
    x[c(seq(1, n, 2), seq(2, n, 2)), ] <- x
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1, ] <- D.le %*% x[1:2, ]
    tmp[n, ] <- D.re %*% x[(n - 1):n, ]
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
