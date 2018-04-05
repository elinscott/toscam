"drape.plot" <- function (x, y, z, z2 = NULL, col = tim.colors(64), zlim = range(z, na.rm=TRUE), zlim2 = NULL, add.legend = TRUE, horizontal = TRUE, theta = 30, phi = 20, ...) 
{
#
# Thanks to JiHO for making corrections and useful extensions to this function
#
# if x is a list, discard y and z and extract them from x
	if (is.list(x)) {
		z <- x$z
		y <- x$y
		x <- x$x
	}
	NC <- length(col)
	M <- nrow(z)
	N <- ncol(z)

	
# if z2 is passed ( values for coloring facets ) use it
# if not use the z matrix that is also used to draw the 
# perspective plot. 
	if (!is.null(z2)) {
		M2 <- nrow(z2)
		N2 <- ncol(z2)
		if ((M != M2) | (N != N2)) {
			stop("draping matrix dimensions must match z")
		}
	}
	else {
		z2 <- z
	}

# if zlim2 has not been passed, set reasonable limits.
# if z2 is passed, set it to the range of z2
# if z2 is not passed, z2=z so we set it to the range of z (equal to zlim)
	if (is.null(zlim2)) {
		 zlim2 <- range(c(z2), na.rm = TRUE)
	}

# determine the colors for facets based on z2, the color scale and 
# the zlim2 z  limits
	zcol <- drape.color(z2, col = col, zlim = zlim2)

# draw filled wireframe and save perspective information
	pm <- persp(x, y, z, theta = theta, phi = phi, col = zcol, zlim=zlim, ...)

# Well, any guesses what this does?
# Note that zlim2 defines limits of color scale

	if (add.legend) {
		image.plot(zlim = zlim2, legend.only = TRUE, col = col, 
			horizontal = horizontal)
	}

# return pm if an assignment is made (see help file)
	invisible(pm)
}

