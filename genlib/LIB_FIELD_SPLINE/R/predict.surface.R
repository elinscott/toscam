"predict.surface" <-
function (object, grid.list = NA, extrap = FALSE, chull.mask=NA, 
     nx = 80, ny = 80, xy=c(1,2),order.variables="xy",...) 
{

#
# without grid.list 
# default is 80X80 grid on first two variables 
# rest are set to median value of x. 
#

  if (is.na(grid.list)[1]) {
        if (is.null(object$x)){ 
          stop("Need a an X matrix in the output object")}
        grid.list<- fields.x.to.grid(object$x, nx=nx, ny=ny,xy=xy)
  }

#
# create grid     
    xg<- make.surface.grid( grid.list)

# make predictions ...

    z<- predict(object, xg, ...)
#
# coerce to the plotting format ($x $y $z etc.) gridding info
# is in grid.list

    out <-  as.surface( grid.list,z, order.variables= order.variables)

# 
# if extrapolate is FALSE set all values outside convex hull to NA
#
    if (!extrap) {
        if (is.na(chull.mask) ) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE) == 
            0] <- NA
    }
#

    out
}

