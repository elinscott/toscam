Krig.null.function<- function( x, Z=NULL,drop.Z=FALSE, m){

# default function to create matrix for fixed part of model
#  x, Z, and drop.Z are required
#  Note that the degree of the polynomial is by convention (m-1)
# returned matrix must have the columns from Z last!
# 
  if( is.null( Z)| drop.Z){
    
     return(fields.mkpoly( x, m=m))}
  else{
     return(cbind(fields.mkpoly( x, m=m),Z)) }

}
