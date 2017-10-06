"Rad.simple.cov" <-
function (x1, x2, p = 1, with.log = TRUE, with.constant = TRUE, 
C= NA,marginal=FALSE) {

    if( marginal) {
     return( rep( 1, nrow( x1)) )}

    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)

    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    m <- (d + p)/2
    
    temp <- rdist( x1,x2)

    if (with.constant) {
        Amd <- radbas.constant(m, d)
    }
    else {
        Amd <- 1
    }

    if ((d%%2 == 0) & (with.log)) {
        temp <- Amd* ifelse(temp < 1e-10, 0, temp^(p/2) * log(temp))}
    else{
         temp <- Amd*temp^(p)}

#
# 
    if( is.na(C[1]) ){ 
     return( temp)}
    else{
     return( temp%*%C) } 


}
