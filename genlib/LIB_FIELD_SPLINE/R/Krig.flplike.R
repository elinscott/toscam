"Krig.flplike" <-
function (lam, obj) 
{

#  - log profile likelihood for lambda
# See section 3.4 from Nychka  Spatial Processes as Smoothers paper. 
# for equation and derivation

    lD <- obj$matrices$D * lam
    nn<- length( lD)
# 
    num<-  log( sum( (obj$matrices$u**2) * lD/(1 + lD) ) )/(nn-obj$nt)

# note subtle differences between den and RSS in Krig.fgcv 

#   log determinant of I-A restricted to nonzero eigenvalues
#   log det is the sum of the logs of the eigenvalues

    den<- sum( ifelse( lD>0, log(lD/(1 + lD)),0))
# NOTE: minus the likelihood! 
   .5* ( num-den) 
 
}
