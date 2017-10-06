"plot.vgram.matrix" <-
function(x,...){
# check if just radil distance has been used for vgram
collapse<- is.null( x$ind)
#
if( !collapse){

nx<- max( x$ind[,1])
ny<- max(  x$ind[,2])
temp<-  matrix( NA,nrow=nx+1, ncol=ny+1)
temp[ x$ind+1] <- x$vgram
image.plot( 0:nx, 0:ny, temp, xlab="X", ylab="Y",...)
}
else( plot( x$d, x$vgram))

}
