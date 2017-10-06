"interp.surface.grid" <-
function( obj, grid.list){
x<- grid.list$x
y<- grid.list$y
M<- length(x)
N<- length( y)
out<- matrix( NA, nrow=M, ncol=N)
for(  i in 1:M){
out[i,] <- interp.surface( obj, cbind( rep( x[i],N),y))
}
list( x=x, y=y, z=out)
}

