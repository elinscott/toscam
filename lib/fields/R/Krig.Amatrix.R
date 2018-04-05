
Krig.Amatrix<- function (object, x0 = object$x, lambda=NULL)
{

  if( is.null( lambda)){ lambda<- object$lambda}

  M<- nrow(object$xM)
  N<- nrow( x0)

# create output matrix 
  out<- matrix( NA,N,M)

#
# loop through unique data locations predicting response 
# using unit vector
# NOTE that the y vector has already been collapsed onto means. 
#

    for( k in 1: M){
     ytemp<- rep( 0,M)
     ytemp[k] <- 1
     out[,k] <- predict(object, x= x0, yM= ytemp, lambda=lambda)
    } 

    return(out)
}
 

