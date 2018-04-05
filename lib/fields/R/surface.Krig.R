"surface.Krig" <-
function(obj, grid.list = NA, extrap = FALSE, graphics.reset = NULL, xlab = NULL,
	ylab = NULL, main = NULL, zlab = NULL, zlim = NULL, levels = NULL,
	type = "C", nx=80, ny=80, ...)
{
	## modified so that you can give main, and ylab as arguments
	## in ... and have them passed correctly
	out.p <- predict.surface(obj, 
             grid.list = grid.list, extrap = extrap, 
               nx=nx,ny=ny, drop.Z=TRUE)
	if(!is.null(ylab))
		out.p$ylab <- ylab
	if(!is.null(xlab))
		out.p$xlab <- xlab
	if(!is.null(zlab))
		out.p$zlab <- zlab
	if(!is.null(main))
		out.p$main <- main
	##    else
	##      out.p$main <- NULL
	plot.surface(out.p, type = type, graphics.reset = graphics.reset, 
		levels = levels, zlim = zlim, ...)
	invisible()
}
