"cat.matrix" <-
function (mat, digits = 8) 
{
    nc <- ncol(mat)
    temp <- matrix(match(c(signif(mat, digits)), unique(c(signif(mat, 
        digits)))), ncol = nc)
    temp2 <- format(temp[, 1])
    if (nc > 1) {
        for (k in 2:nc) {
            temp2 <- paste(temp2, temp[, k], sep = "X")
        }
    }
    match(temp2, unique(temp2))
}
