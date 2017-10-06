"quilt.plot" <-
function(x,y,z,nrow=64, ncol=64, grid=NULL, add.legend=TRUE,add=FALSE,...){


x<- as.matrix(x)

  if( ncol(x)==2){
   z<- y}

  if( ncol(x)==1){
   x<- cbind( x,y)}
  if( ncol(x)==3){
  z<- x[,3]
  x<- x[,1:2]
}

# at this point x should be a 2 column matrix of x-y locations
#  z is a vector or one column matrix of the z values.

#discretize data
as.image( z, x=x, nrow=nrow, ncol=ncol, na.rm=TRUE)-> out.p

#plot it 
if( add.legend){
image.plot( out.p,col=tim.colors(64),add=add,...)
}
else{
image(out.p, col=tim.colors(64),add=add,...)
}


}

