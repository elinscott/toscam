Krig.make.W<- function( out, verbose=FALSE){

   if( verbose) { cat( "W", fill=TRUE) ; print( out$W)}
    
   if( out$nondiag.W) {
#
# create W from scratch or grab it from passed object

        if( is.null(out$W)){
         if( verbose){ print( out$wght.function.name)}
         W<- do.call( out$wght.function.name, 
                        c( list( x=out$xM), out$args.wght))
#
# possibly adjust W based on diagional weight terms 
#          
         W<- sqrt( out$weightsM)* t( sqrt(out$weightsM)*W) }
        else{ 
          W<-  out$W}
#
# symmetric square root
      eigen( W, symmetric=TRUE)-> temp
      W2<- temp$vectors%*% diag( sqrt( temp$values))%*% t( temp$vectors)
     return( list( W=W, W2=W2)) }
   
   else{
#
#  If W is diagonal include  as a vector. 
#  Subsequent multiplies use %d*% to handle the diagonal and nondiagonal cases
#  together   e.g.  W%d*%yM will worrk for both diag and nondiag W's. 
#
      return(
             list( W= out$weightsM, W2=sqrt(out$weightsM))
             )}

   
}
