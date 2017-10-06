"world" <-
function (ylim = c(-90, 90), xlim = NULL, add = FALSE, asp = 1, 
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", eps = 0.1, 
    shift = FALSE, ...) 
{
    if (shift) {
ind1<- !is.na(world.dat$x)
ind2<- (world.dat$x<0)
# shift coordinates
        world.dat$x[ind2&ind1] <- world.dat$x[ind2&ind1] + 360
# "pick up pen" at the new edges by inserting NA's
        world.dat$x[(world.dat$x <= eps | world.dat$x >= (360 - 
            eps))&ind1] <- NA
    }
    if (is.null(xlim)) {
        if (shift) {
            xlim <- c(0, 360)
        }
        else {
            xlim <- c(-180, 180)
        }
    }
    if (!add) {
        plot(world.dat, ylim = ylim, xlim = xlim, type = "n", 
            xaxt = xaxt, yaxt = yaxt, xlab = xlab, ylab = ylab, 
            asp = asp, ...)
    }
    lines(world.dat, err = -1, ...)
    invisible()
}

