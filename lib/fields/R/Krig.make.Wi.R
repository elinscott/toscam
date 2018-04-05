Krig.make.Wi<- function( out, verbose=FALSE){
#
# If a weight matrix has been passed use it.
#  
# Note that in either case the weight matrix assumes that 
# replicate observations have been collapses to the means. 
#
   if( out$nondiag.W) {
      eigen( out$W, symmetric=TRUE)-> temp
      Wi<- temp$vectors%*% diag( 1/( temp$values))%*% t( temp$vectors)
      W2i<- temp$vectors%*% diag( 1/sqrt( temp$values))%*% t( temp$vectors)
     return( list( Wi=Wi, W2i=W2i)) }
   
   else{
#
#  These are created only for use with default method to stay 
# consistent with nondiagonal elements. 
      return(
             list( Wi= 1/out$weightsM, W2i= 1/sqrt(out$weightsM) )
             )}

   
}
