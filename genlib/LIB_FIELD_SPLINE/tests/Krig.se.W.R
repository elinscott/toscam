library( fields)
# tests of predict.se using 
# off diag weight matrix for obs (W)

options( echo=FALSE)


# a nasty example with off diagonal weights.

set.seed(123)

N<- 50 
x<- matrix( runif( N*2), N,2)
y<- rnorm( N)*.2 + 2*x[,1]**2 +  x[,2]**2 


weights<- runif(N)*10
x0<- cbind( c(.1,.2,.6,.65,.8), c(.05,.5,.73,.9,.95))



temp.wght<- function(x, alpha=.3){
  Exp.cov( x, theta=.1) }

Krig( x,y, cov.function=Exp.cov,weights=weights,
     wght.function= "temp.wght")-> out
Krig( x,y, cov.function=Exp.cov,weights=weights,W= out$W)-> out2


# direct calculation test for A matrix
#

Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%y, predict( out, x0),tag="Amatrix vs. predict")

# now find se.

W2<-out$W2
W<- out$W

Sigma<- out$rhohat*Exp.cov( out$x,out$x)
temp0<- out$rhohat*(Exp.cov( x0, x0))
S1<- out$rhohat*Exp.cov( out$x, x0)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

Sigma.obs<-  Krig.make.Wi( out)$Wi
Sigma.obs <- Sigma.obs* (out$shat.MLE**2) 

temp1<-  A%*%S1
temp2<- A%*% ( Sigma.obs+ Sigma)%*% t(A)
look<- temp0 - t(temp1) - temp1 +  temp2


#compare to 
# diagonal elements

test<- predict.se( out, x= x0) 
test.for.zero( sqrt(diag(  look)), test,tag="Marginal predict.se")


test<- predict.se( out, x= x0, cov=TRUE)
test2<- predict.se( out2, x= x0, cov=TRUE)
test.for.zero( look, test,tag="Full covariance predict.se")
test.for.zero( look, test2,tag="Full covariance predict.se explicit W")

cat( "all done", fill=TRUE)
options( echo=TRUE)
