"Wimage.cov" <-
function( Y=NULL,ind=NULL, H.obj, find.row=FALSE ){

# local multipcation by H function
   Mult.H<-function( u, H.obj){
            IND<- H.obj$ind0
            u[IND] <- H.obj$H0%*%c(u[IND])
            u*H.obj$H1
         }
#
# This is a shortcut option to find a single row of the covariance matrix
#
    if( find.row){ Y<-1}
#
# if ind is  passed fill complete image with the Y subset       
# 
   if(is.null(ind)){
      tmp<- Y}
    else{
      tmp<- matrix(0,nrow=H.obj$m,ncol=H.obj$n)
      tmp[ind]<- Y
    }
#
# 
     hold<-  Wtransform.image( tmp, cut.min=H.obj$cut.min, inv=TRUE,transpose=TRUE)
#
     hold<- Mult.H( hold, H.obj=H.obj) # Note if H not symmetric this should be 
                                        # multiplication by transpose.
     hold<- Mult.H( hold, H.obj=H.obj)
     Wtransform.image( hold, cut.min=H.obj$cut.min, inv=TRUE)
}
