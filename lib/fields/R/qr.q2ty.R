"qr.q2ty" <-
function (qr, y) 
{
    if (!is.matrix(y)) {
        y <- as.matrix(y)
    }
    dy <- dim(y)
    dq <- dim(qr$qr)
    rank <- qr$rank
    if (dy[1] != dq[1]) 
        stop("y and qr$qr should have same number of rows")
    qr.qty(qr, y)[(rank + 1):dy[1], ]
}
