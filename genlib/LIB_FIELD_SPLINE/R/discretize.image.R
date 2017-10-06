"discretize.image" <-
function(x, m = 64, n = 64,  grid = NULL, 
         expand = c(1, 1), boundary.grid=FALSE)
{
	#
	# set up discretized grid based on x
	#

	out <- list()

	if(length(expand) == 1)
		expand <- rep(expand, 2)

	if(is.null(grid)) {
               grid<- list()
               xr <- range(x[, 1])
	 	deltemp <- (xr[2] - xr[1]) * (expand[1] - 1) * 0.5
		grid$x <- seq(xr[1] - deltemp, xr[2] + deltemp,  , m)

		yr <- range(x[, 2])
		deltemp <- (yr[2] - yr[1]) * (expand[2] - 1) * 0.5
		grid$y <- seq(yr[1] - deltemp, yr[2] + deltemp,  , n)}

        

# find cut points for boundaries assuming midpoints
       if( !boundary.grid){        
	xcut <- fields.convert.grid(grid$x)
	ycut <- fields.convert.grid(grid$y)}

# cut points given boundaries
       else{
          xcut<- grid$x
          ycut<- grid$y}

# locate bin ids for each location
        
	out$index <- cbind(cut(x[, 1], xcut), cut(x[, 2], ycut))

	out$m <- length(xcut)-1
	out$n <- length(ycut)-1

	out$grid <- grid
        if( !boundary.grid ){
#discretized locations          
	out$loc <- cbind( grid$x[out$index[, 1]], grid$y[out$index[, 2]] )}
        else{
          out$loc<- NA}

        out
}
