"interp.surface" <-
function (obj, loc)
{

# obj is a surface object like the list for contour or image.
# loc a matrix of 2 d locations 

# Thanks to Steve Koehler for this 
# compact version of bilinear interpolation on a 2-d grid

x<- obj$x
y<- obj$y
x.new<- loc[,1]
y.new<- loc[,2]
z<- obj$z
nx<- length( x)
ny<- length( y)
# x and y are the grid values assumed to be equally spaced
# z is a matrix of grid values 

# rescale the points to be interpolated with respect to
# the grid spacing. 

   lx <- ((length(x) - 1) * (x.new - min(x)))/diff(range(x)) + 1
   ly <- ((length(y) - 1) * (y.new - min(y)))/diff(range(y)) + 1
#
# nearest grid point that is less than the new x and y values 
# set to NA any points that are outside of grid (but see further code below)

   lx1 <- floor (lx); lx1[lx1 < 1 | lx1 >= nx] <- NA
   ly1 <- floor (ly); ly1[ly1 < 1 | ly1 >= ny] <- NA

#
# difference ( in scale of grid spacing) between the
#  new x, y values and the lower grid values 

   ex <- lx - lx1
   ey <- ly - ly1

# bilinear interpolation finds simple weights based on the
# the four corners of the grid box containing the new 
# points. 

# fix up weights to handle the case when loc are equal to 
# last grid point.  These have been set to NA above.

   max.x<- max(x)
   max.y<- max( y) 

   ex[ x.new == max.x] <- 1.0
   ey[ y.new == max.y] <- 1.0
   lx1[ x.new == max.x] <- nx-1
   ly1[ y.new == max.y] <- ny-1

   return (

      z[cbind(lx1,     ly1)]     * (1 - ex) * (1 - ey) +
      z[cbind(lx1 + 1, ly1)]     * ex       * (1 - ey) +
      z[cbind(lx1,     ly1 + 1)] * (1 - ex) * ey       +
      z[cbind(lx1 + 1, ly1 + 1)] * ex       * ey           )
}

