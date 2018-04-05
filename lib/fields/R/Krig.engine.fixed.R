"Krig.engine.fixed" <-
function(out, verbose=FALSE, lambda=NA){
# 
# Model:
#     Y_k=  f_k + e_k
#  var( e_k) = sigma^2/W_k
#
#   f= Td + h 
#    T is often a low order polynomial 
#   E(h)=0    cov( h)= rho *K
#
# let M = (lambda W^{-1} + K)
# the surface estimate  depends on coefficient vectors d and c
#    The implementation in Krig/fields is that K are the 
#    cross covariances among the observation locations and the knot locations
#    H is the covariance among the knot locations. 
#    Thus if knot locs == obs locs we have the obvious collapse to 
#    the simpler form for M above.  
#
#   With M in hand ...
#
#   set 
#   d =  [(T)^t M^{-1} (T)]^{-1} (T)^t M^{-1} Y
#  this is just the generalized LS estimate for d 
#
#   lambda= sigma**2/rho 
#  the estimate for c is  
#   c=  M^{-1}(y - Td)
#
# This particular numerical strategy takes advantage of 
# fast Cholesky factorizations for positive definite matrices
# and also provides a seamless framework for sparse matrix implementations
#

  if( is.na( lambda)) lambda<- out$lambda

  call.name<- out$cov.function.name

  if( !out$knot.model){
####################################################
# case of knot locs == obs locs  out$knots == out$xM
####################################################

# create T matrix
  Tmatrix <-  do.call(out$null.function.name, 
                        c(out$null.args, list(x=out$knots, Z=out$ZM) ) )

  if( verbose){
     cat("Tmatrix:", fill=TRUE)
     print( Tmatrix)}
  np <- nrow( out$knots)
  nt <- ncol( Tmatrix)
# form K 
     tempM<-  do.call(
             call.name,
            c(out$args, list(x1 = out$knots, x2 = out$knots))  )
# form M 

     diag(tempM) <- (lambda/ out$weightsM) + diag(tempM)
#
# find cholesky factor  
#  tempM = t(Mc)%*% Mc
#  V=  Mc^{-T}   

# call cholesky but also add in the args supplied in Krig object. 

  do.call("chol",c(list( x=tempM), out$chol.args))-> Mc

 VT<- forwardsolve( Mc, x=Tmatrix, transpose=TRUE, upper.tri=TRUE)
 qr( VT)-> qr.VT

#
# now do generalized least squares for d 
# and then find c. 

 d.coef<- qr.coef( qr.VT, forwardsolve( Mc, transpose=TRUE, out$yM, upper.tri=TRUE) )

 if( verbose){
  print( d.coef)}

 c.coef<- forwardsolve( Mc, transpose=TRUE, out$yM- Tmatrix%*% d.coef, upper.tri=TRUE)
 c.coef<- backsolve( Mc,c.coef)

# return all the goodies,  include lambda as a check because 
# results are meaningless for other values of lambda 

return(
  list(  qr.VT = qr.VT, d=c(d.coef), c=c(c.coef), 
         Mc = Mc, decomp = "cholesky", 
            nt=nt, np=np,lambda.fixed=lambda)
 )}
 else { 
####################################################
# case of knot locs != obs locs  
####################################################
# create weighted T matrix

  Tmatrix <-  do.call(out$null.function.name, 
                        c(out$null.args, list(x=out$xM, Z= out$ZM) ) )

  nt <- ncol( Tmatrix)
  np <- nrow( out$knots) +nt
# form H
   H <-  do.call(
            call.name,
            c(out$args, list(x1 = out$knots, x2 = out$knots) ) )

# form K matrix 
    K<-  do.call(
            call.name,
            c(out$args, list(x1 = out$xM, x2 = out$knots) ) )
# 
  
     Mc<- do.call("chol",
       c(list( x=t( K) %*%(out$weightsM * K) + lambda*H ), out$chol.args) )

# weighted Y
    wY<-  out$weightsM* out$yM

    temp0<- t(K)%*% (out$weightsM * Tmatrix)
    temp1<- forwardsolve( Mc, temp0, transpose=TRUE, upper.tri=TRUE)

    qr.Treg<- qr(t(Tmatrix)%*%(out$weightsM*Tmatrix) - t(temp1)%*%temp1)
    
    temp0<- t(K)%*% wY
    temp3<- t(Tmatrix)%*% wY - 
                 t(temp1)%*%forwardsolve( Mc, temp0, transpose=TRUE , upper.tri=TRUE)

     d.coef<- qr.coef( qr.Treg, temp3)

     temp1<- t(K)%*% (wY - out$weightsM * (Tmatrix)%*% d.coef)
     c.coef<- forwardsolve( Mc, transpose=TRUE, temp1, upper.tri=TRUE)
     c.coef<- backsolve( Mc, c.coef)

  list(  qr.Treg=qr.Treg, d=c(d.coef), c=c(c.coef), 
         Mc= Mc, decomp = "cholesky.knots", 
            nt=nt, np=np,lambda.fixed=lambda)

}
#
# should not get here. 
#
}

