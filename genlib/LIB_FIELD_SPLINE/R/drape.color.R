"drape.color" <- function (z, col=tim.colors(64),zlim=NULL,transparent.color="white") 
{
# slight increase for range 
	eps <- 1e-7
# range if zlim not supplied
	if( is.null(zlim)){
		zlim <- range( c(z), na.rm=TRUE)
	} 
# set any values outside of range to NA ( i.e. the transparent.color)
	z[ (z< zlim [1]) | (z> zlim [2]) ] <- NA
	
	NC <- length( col)
	M <- nrow(z)
	N <- ncol(z)
# find average z value for a facet
	zm <- ( z[1:(M-1),1:(N-1)] + z[2:M, 1:(N-1)] + z[1:(M-1),2:N] + z[2:M,2:N] )/4
# spacing for grid to assign  colors
# 1+eps included so that if z== zlim[2] it gets a color
	dz <- (zlim[2]*(1+eps)- zlim[1])/NC
# discretize the z-value 
	zcol <- floor( (zm -zlim[1])/dz  +1)
	
# assign colors if in range 1:NC otherwise use the transparent color
	ifelse(zcol > NC, transparent.color, col[zcol])

}

