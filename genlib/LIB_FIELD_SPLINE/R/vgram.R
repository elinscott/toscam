"vgram" <-
function(loc, y, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL, 
N = NULL, breaks = NULL)
{
	if(is.null(id)) {
		n <- nrow(loc)
		ind <- rep(1.:n, n) > rep(1.:n, rep(n, n))
		id <- cbind(rep(1.:n, n), rep(1.:n, rep(n, n)))[ind,  ]
	}
	##
	## OK here it is :
	##
	if(is.null(d)) {
		loc <- as.matrix(loc)
		if(lon.lat) {
			d <- rdist.earth(loc)[id]
		}
		else {
			d <- rdist(loc, loc)[id]
		}
	}
	vg <- 0.5 * (y[id[, 1.]] - y[id[, 2.]])^2.
	call <- match.call()
	##
	##
	##
	if(is.null(dmax)) {
		dmax <- max(d)
	}
	od <- order(d)
	d <- d[od]
	## Replace next two statements
	## vgram <- vgram[od]
	## ind <- d <= dmax & !is.na(vgram)
	vg <- vg[od]
	ind <- d <= dmax & !is.na(vg)
	##
	## binned  variogram if breaks supplied
	##
	out <- list(d = d[ind], vgram = vg[ind], call = call)
	if(!is.null(breaks) | !is.null(N)) {
		out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
	}
	out
}
