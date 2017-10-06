"Krig.engine.default" <-
function(out, verbose=FALSE){
# 
# matrix decompositions for computing estimate

#
# Computational outline ( "." is used for subscript)
#
# The form of the estimate is
#    fhat(x) = sum phi.j(x) d.j  + sum psi.k c.k
#
# the {phi.j} are the fixed part of the model usually low order polynomials
# and is also referred to as spatial drift.
#
# the {psi.k} are the covariance function evaluated at the unique observation
# locations.  if xM.k is the kth unique location psi.k(x)= k(x, xM.k)
# xM is also out$knots in the code below.
#
# the goal is find decompositions that facilitate rapid solution for
# the vectors d and c. 
#

#  First outline calculations with equal weights
#  T the fixed effects regression matrix  T.ij = phi.j(xM.i)
#  K the covariance matrix for the unique locations

# From the spline literature the solution solves the well known system 
# of two eqautions:

#    -K( yM - Td - Kc) + lambda *Kc = 0
#     -T^t ( yM-Td -Kc) =0

# Divide through by K and substitute, these are equivalent to 

#  -1-      -( yM- Td - Kc) + lambda c = 0
#  -2-      T^t c = 0       



#
#  A QR decomposition is done for   T= (Q.1,Q.2)R
#   by definition  Q.2^T T =0
  
#  equation  -2- can be thought of as a constraint
# with  c= Q.2 beta2
# substitute in  -1-  and multiply through by Q.2^T

#      -Q.2^T yM  + Q.2^T K Q.2 beta2  + lambda beta2 = 0
#
#   Solving
#   beta2 = {Q.2^T K Q.2 beta2  + lambda I )^ {-1} Q.2^T yM
# 
#  eigenvalues and eigenvectors found for M= Q.2^T K Q.2 
#
#
# From -1-  Td = yM - Kc - lambda c
#      (Q.1^T) Td =  (Q.1^T) ( yM- Kc)    
#
#   ( lambda c is zero by constraint)
#   
#   so Rd = (Q.1^T) ( yM- Kc)
# use qr functions to find d. 
# 


        Tmatrix<- do.call(out$null.function.name, 
                           c(out$null.args, list(x=out$xM,Z=out$ZM))  )
        if( verbose){
             cat(" Model Matrix: spatial drift and Z", fill=TRUE)
             print( Tmatrix)
        }

# Tmatrix premultiplied by sqrt of wieghts
        Tmatrix<- out$W2%d*%Tmatrix 

        qr.T <- qr(Tmatrix )

#
#verbose block
        if (verbose) {
            cat( "first 5 rows of qr.T$qr",fill=TRUE)
            print(qr.T$qr[1:5,])
        }
#
# find  Q_2 K Q_2^T  where K is the covariance matrix at the knot points
#
        
        tempM<-   t(
           out$W2 %d*% 
           do.call(out$cov.function.name, 
               c(out$args, list(x1 = out$knots, x2 = out$knots)))
                    )

        tempM <- out$W2 %d*% tempM
        
        tempM <- qr.yq2(qr.T, tempM)
        tempM <- qr.q2ty(qr.T, tempM)

    np <- nrow(out$knots)
    nt <- (qr.T$rank)

    if (verbose) {
        cat("np, nt", np, nt, fill = TRUE)}

#
# Full set of decompositions for 
# estimator for nonzero lambda

            tempM <- eigen(tempM, symmetric=TRUE)
            D <- c(rep(0, nt), 1/tempM$values)
#
# verbose block
            if (verbose) {
                cat("eigen values:", fill = TRUE)
                print(D)
            }
#
# Form the matrix decompositions and transformed data vector used to 
# evaluate the solution, GCV, REML  at different lambdas
#
            G <- matrix(0, ncol = np, nrow = np)
        
            G[(nt + 1):np, (nt + 1):np] <- tempM$vectors

            G <- G * matrix(D, ncol = np, nrow = np, byrow = TRUE)

            u <- c( rep(0, nt),
                t(tempM$vectors) %*% 
                    qr.q2ty(qr.T, c(out$W2%d*%out$yM ) ))
#
# verbose block
            if (verbose){
                print(u)
                print(out$pure.ss)
            }
#
            return( list(u = u, D = D, G = G, qr.T = qr.T, 
                decomp = "WBW", V = tempM$vectors, nt=nt, np=np))
}

