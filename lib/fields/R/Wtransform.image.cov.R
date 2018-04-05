"Wtransform.image.cov" <-
function (ind1, ind2=ind1,Y, cov.obj) 
{

#
# define useful local function that does the H multiplcation
# taking advantage of the block diagonal form 
#  

    Mult.H <- function(u, cov.obj) {
        IND <- H.obj$ind0
        u[IND] <- H.obj$H0 %*% c(u[IND])
        u * H.obj$H1
    }

         IND <- cov.obj$ind0 
         CUT<- cov.obj$cut.min         

      if( missing(ind1)){ 
      # do multiplcation for full matrix. 
         if( (nrow(Y)!=cov.obj$m) | (ncol(Y)!=cov.obj$m)){ 
             stop("bad dimensions for Y")}
         hold <- Wtransform.image(Y,cut.min=CUT,inv=TRUE,transpose = TRUE)
         hold <- Mult.H(hold, cov.obj)
         hold <- Mult.H(hold, cov.obj)
         return(
              Wtransform.image(hold, cut.min = CUT, inv = TRUE))
        }

        else{
        # multiplication for subset of indices

         temp <- matrix(0, nrow = cov.obj$m, ncol = cov.obj$n)
         temp[ind2] <- Y 
#
         hold <- Wtransform.image(temp,cut.min=CUT,inv = TRUE,transpose=TRUE)
         hold <- Mult.H(hold, cov.obj)
         hold <- Mult.H(hold, cov.obj)
         return( 
              Wtransform.image(hold, cut.min = CUT,inv =TRUE)[ind1] )
         }

}

