
library( fields)

# tests of predict.se
# against direct linear algebra 

options( echo=FALSE)

x0<- expand.grid( c(-8,-4,0,20,30), c(10,8,4,0))


Krig( ozone$x, ozone$y, cov.function = "Exp.cov", theta=50)-> out


# direct calculation
Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%ozone$y, predict( out, x0),tag="Amatrix vs. predict")

Sigma<- out$rhohat*Exp.cov( ozone$x, ozone$x, theta=50)
S0<- out$rhohat*c(Exp.cov( x0, x0, theta=50))
S1<- out$rhohat*Exp.cov( out$x, x0, theta=50)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

look<- S0 - t(S1)%*% t(A) - A%*%S1 +  
       A%*% ( Sigma + diag(out$shat.MLE**2/out$weightsM))%*% t(A)
#
#compare to 
# diagonal elements


test2<- predict.se( out, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predict.se")

out2<- Krig( ozone$x, ozone$y, cov.function = "Exp.cov", theta=50,
            lambda=out$lambda)

test2<- predict.se( out2, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predict.se fixed ")

test<- predict.se( out, x= x0, cov=TRUE)
test.for.zero( look, test,tag="Full covariance predict.se")


# simulation based.

set.seed( 333)

sim.Krig.standard( out, x0,M=4e3)-> test

var(test)-> look

predict.se( out, x=x0)-> test2
mean( diag( look)/ test2**2)-> look2
test.for.zero(look2, 1.0, tol=1e-2, tag="Marginal standard Cond. Sim.")

predict.se( out, x=x0, cov=TRUE)-> test2

# multiply simulated values by inverse square root of covariance
# to make them white

eigen( test2, symmetric=TRUE)-> hold
hold$vectors%*% diag( 1/sqrt( hold$values))%*% t( hold$vectors)-> hold
cor(test%*% hold)-> hold2
# off diagonal elements of correlations -- expceted values are zero. 

abs(hold2[ col(hold2)> row( hold2)])-> hold3

test.for.zero(   mean(hold3), 0, relative=FALSE, tol=.02,
          tag="Full covariance standard Cond. Sim.")


# test of sim.Krig.grid.R
#
# first create and check a gridded test case. 


data( ozone2)
as.image(ozone2$y[16,], x= ozone2$lon.lat, ncol=24, nrow=20, 
          na.rm=TRUE)-> dtemp
#
# A useful disctrtized version of ozone2 data
 
x<- cbind(dtemp$x[dtemp$ind[,1]], dtemp$y[dtemp$ind[,2]])
y<- dtemp$z[ dtemp$ind]
weights<- dtemp$weights[ dtemp$ind]

Krig( x, y, Covariance="Matern", 
   theta=1.0, smoothness=1.0, weights=weights) -> out



  set.seed(234)
  ind0<- cbind( sample( 1:20, 5), sample( 1:24, 5))

  x0<- cbind( dtemp$x[ind0[,1]], dtemp$y[ind0[,2]]) 

# an  inline check plot(out$x, cex=2); points( x0, col="red", pch="+",cex=2)

# direct calculation as backup ( also checks weighted case)

Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%out$yM, predict( out, x0),tag="Amatrix vs. predict")

Sigma<- out$rhohat*stationary.cov( 
out$xM, out$xM, theta=1.0,smoothness=1.0, Covariance="Matern")

S0<- out$rhohat*stationary.cov( 
x0, x0, theta=1.0,smoothness=1.0, Covariance="Matern")

S1<- out$rhohat*stationary.cov(
out$xM, x0, theta=1.0,smoothness=1.0, Covariance="Matern")



#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)
 
look<- S0 - t(S1)%*% t(A) - A%*%S1 +
       A%*% ( Sigma + diag(out$shat.MLE**2/out$weightsM) )%*% t(A)

test<- predict.se( out, x0, cov=TRUE)

test.for.zero( c( look), c( test), tag="Weighted case and exact for ozone2 full 
cov", tol=1e-8)


#cat("Conditional simulation test -- this takes some time", fill=TRUE)

# the grid ...

glist<- list( x= dtemp$x, y=dtemp$y)

set.seed( 233)
sim.Krig.grid( out, grid= glist, M=400)-> look

predict.surface.se( out, grid=glist)-> test

look2<- matrix( NA, 20,24)

for(  k in 1:24){
for ( j in 1:20){
look2[j,k] <- sqrt(var( look$z[j,k,], na.rm=TRUE))
}
}


test.for.zero(  1-mean(c(look2/test$z), na.rm=TRUE), 0, relative=FALSE, 
tol=.005, tag="Conditional simulation marginal se for grid")

cat("all done testing predict.se ", fill=TRUE)
options( echo=TRUE)
