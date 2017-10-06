"sim.Krig.standard" <-
function( object, xp, M=1, verbose=FALSE, sigma2=NA, rho=NA) {

# figure out what sigma and rho should be
 
 if( is.na(sigma2)){
   sigma2<- object$best.model[2]}
 if( is.na( rho)){
    rho<- object$best.model[3]} 

#
# check for unique rows of xp
 
  if( any( duplicated(xp))){ 
     stop(" predictions locations should be unique")}

#
# set up various sizes of arrays
     m<- nrow( xp)
     n<- nrow(object$xM) # NOTE: xM _unique_ locations of data
     N<- length( object$y)
     if( verbose) { 
       cat( " m,n,N, sigma2, rho", m, n, N, sigma2, rho,fill=TRUE)}

#transform the new points
  xc <- object$transform$x.center
  xs <- object$transform$x.scale
  xpM <- scale(xp, xc, xs)

# complete set of points for prediction.
# check for replicates and adjust

     x<- rbind( object$xM, xpM)
#
# find indices of all rows of xp that correspond to rows of 
# xM and then collapse x to unique rows. 

     rep.x.info<- fields.duplicated.matrix(x)
     x<- x[ !duplicated(rep.x.info), ]    
     N.full<- nrow( x)

# these give locations in x matrix to reconstruct xp matrix

     xp.ind<- rep.x.info[(1:m)+n]

     if( verbose){ 
       print( N.full) 
       print( x)
     }

     if( verbose){
      cat( "reconstruction of xp from collapsed locations", 
            fill=TRUE)
          print(x[xp.ind,]) }
#
# Sigma is full covariance at the data locations and at prediction points. 
#

    Sigma <-  rho * do.call(
                  object$cov.function.name, 
                   c( object$args, list(x1 = x, x2 = x )))
#
# square root of Sigma for simulating field
# Cholesky is fast but not very stable. 
#
# the following code line is similar to chol(Sigma)-> Scol
# but adds possible additional arguments controlling the Cholesky 
# from the Krig object. 
#

   do.call("chol",c(list( x=Sigma), object$chol.args) ) ->  Schol

#
# output matrix to hold results
 
    N.full<- nrow( x)
    out<- matrix( NA, ncol=m, nrow= M)

#
# find conditional mean field from initial fit
# don't multiply by sd or add mean if this is 
# a correlation model fit. 
# (these are added at the predict step). 

    h.hat<- predict( object, xp ) 

# marginal standard deviation of field. 
    temp.sd<- 1
#
# 
# this is not 1 if Krig  object is a corelation model. 

    if (object$correlation.model) {
      if( !is.na( object$sd.obj[1]) ) {temp.sd<-predict( object$sd.obj, x)}
    }
#
#  Define W2i for simulating errors. 
#

  W2i <- Krig.make.Wi( object)$W2i

  for(  k in 1: M){

# simulate full field
       t(Schol)%*%rnorm(N.full)-> h

# value of simulated field at observations
#
# NOTE: fixed part of model (null space) need not be simulated 
# because the  estimator is unbiased for this part. 
# the variability is still captured because the fixed part
# is still estimated as part of the predict step below

     h.data <- h[1:n]

# expand the values according to the replicate pattern
  
     h.data<- h.data[object$rep.info]

# create synthetic data

     y.synthetic<- h.data + sqrt(sigma2)* W2i%d*%rnorm(N)

# predict at xp using these data
# and subtract from "true" value 
# note that true values of field have to be expanded in the 
# case of common locations between xM and xp. 

     h.true<- (h[xp.ind])
     temp.error<- predict( object, xp, y=y.synthetic, 
                 eval.correlation.model=FALSE) - h.true 

# add the error to the actual estimate  (conditional mean)
# and adjust by marginal standard deviation
     
     out[k,]<-  h.hat + temp.error*temp.sd
}

out
}

