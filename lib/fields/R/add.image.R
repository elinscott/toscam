"add.image" <-
function(xpos,ypos,z,adj.x=.5,adj.y=.5,image.width=.15, 
image.height= NULL, col=tim.colors(256), ...) 
{
   m<- nrow( z)
   n<- ncol(z)

   ucord <- par()$usr
   pin <- par()$pin

# if height is missing scale according to width assuming pixels are 
# square. 

   if( is.null( image.height)){
      image.height<- (n/m)*image.width}

# find grid spacing in user coordinates. 
   dy<- image.width*(ucord[4] - ucord[3])
   dx<- image.height* pin[2]*(ucord[2] - ucord[1])/(pin[1])


#
# dx and dy should have the correct ratio given different different scales 
# and also different aspects to the plot window 
#

# find grid to put image in right place. 
    xs<- seq(0, dx,, m+1) +xpos - adj.x*dx
    ys<- seq(0, dy,,n+1) +ypos - adj.y*dy 

  image( xs, ys, z,add=TRUE, col=col,...)
}

