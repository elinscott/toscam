"in.poly" <-
function (xd, xp, convex.hull = FALSE, inflation=1e-7) 
{
    if (convex.hull) {
        xp <- xp[chull(xp), ]
    }
    nd <- nrow(xd)
    np <- as.integer(nrow(xp))
#
# inflate convex hull slightly to include any points actually on the hull
#
     if( convex.hull){
         xm<- matrix(c(mean(xp[,1]), mean( xp[,2])), nrow=np, ncol=2,byrow=TRUE)
         xp<- (xp- xm)*( 1+ inflation) + xm
}

    .Fortran("inpoly", nd = as.integer(nd), as.single(xd[, 1]), 
        as.single(xd[, 2]), np = np, as.single(xp[, 1]), as.single(xp[, 
            2]), ind = as.integer(rep(-1, nd)), PACKAGE="fields")$ind
}
