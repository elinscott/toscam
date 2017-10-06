"Wimage.sim" <-
function( H.obj){

# local function for (fast or sparse) multiplication by H 
   Mult.H<-function( u, H.obj){
            IND<- H.obj$ind0
            u[IND] <- H.obj$H0%*%c(u[IND])
            u*H.obj$H1
         }
# 
      tmp<- matrix(rnorm( H.obj$m*H.obj$n),nrow=H.obj$m,ncol=H.obj$n)
#
# 
     tmp<-  Wtransform.image( tmp, cut.min=H.obj$cut.min, inv=TRUE,transpose=TRUE)
     Wtransform.image( tmp, cut.min=H.obj$cut.min, inv=TRUE)
}
