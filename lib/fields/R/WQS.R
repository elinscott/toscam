"WQS" <-
function (x) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(-1, 3, 3, -1)
    D.le <- matrix(c(2, 3, -1, 2, -3, 1), nrow = 2, ncol = 3, 
        byrow = TRUE)
    D.re <- matrix(c(-1, 3, 2, -1, 3, -2), nrow = 2, ncol = 3, 
        byrow = TRUE)
    D.rough <- c(-1, 3, -3, 1)
    n <- dim(x)[1]
    m <- dim(x)[2]
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:2, ] <- D.le %*% x[1:3, ]
    tmp[(n - 1):n, ] <- D.re %*% x[(n - 2):n, ]
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
