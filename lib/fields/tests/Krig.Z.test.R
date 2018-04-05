library(fields)
#
#
#  test of fixed lambda case
#  Check against linear algebra
#

options( echo=FALSE)

#cat("A very nasty case with knots and weights",fill=TRUE)

set.seed(123)
x<- matrix( runif( 30), 15,2)
Z<- matrix( rnorm(30), 15,2)
y<- rnorm( 15)*.01 + 5*(x[,1]**3 +  (x[,2]-.5)**2) +  (Z[,1] +Z[,2])*.001
knots<- x[1:5,]
#weights<- runif(15)*10

# first without knots compare default to fixed

Krig( x,y,Z=Z, cov.function=Exp.cov, give.warnings=FALSE)-> out.new

Krig( x,y,Z=Z, cov.function=Exp.cov,lambda=1)-> out.new2


##########
## compute test using linear algebra

K<- Exp.cov( x,x)
lambda<-1
M<- (lambda* diag(nrow( x)) + K)
T<- cbind( rep(1,15), x, Z)
temp.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*% y)
temp.c<- solve( M)%*% ( y - T%*% temp.d)

# test for d coefficients
test.for.zero( out.new2$d, temp.d, tag=" d coef")
# test for c coefficents
test.for.zero( out.new2$c, temp.c, tag="c coef" )


####### testing predict function 
hold2<- predict( out.new2, x=x, Z=Z, just.fixed=TRUE)
hold3<- predict( out.new2, x=x, Z=Z, drop.Z=TRUE)
hold4<- predict( out.new2, x=x, Z=Z, drop.Z=TRUE, just.fixed=TRUE)

hold<-T%*%temp.d
test.for.zero( hold, hold2, tag="predict for null fixed" )

hold<-T[,1:3]%*%temp.d[1:3] + K %*% temp.c
test.for.zero( hold, hold3, tag="predict for null spatial"  )

hold<-T[,1:3]%*%temp.d[1:3]
test.for.zero( hold, hold4, tag="predict for null drift" )

######tests where coefficients  are recomputed from object
hold2<- predict( out.new,y=y, lambda=1.0, x=x, Z=Z, just.fixed=TRUE)
hold3<- predict( out.new,y=y, lambda=1.0, x=x, Z=Z, drop.Z=TRUE)
hold4<- predict( out.new, y=y, lambda=1.0, x=x, Z=Z, 
                      drop.Z=TRUE, just.fixed=TRUE)

hold<-T%*%temp.d
test.for.zero( hold, hold2, tag="predict for null fixed" )

hold<-T[,1:3]%*%temp.d[1:3] + K %*% temp.c
test.for.zero( hold, hold3, tag="predict for null spatial" )

hold<-T[,1:3]%*%temp.d[1:3]               
test.for.zero( hold, hold4, tag="predict for null drift " )



###knots case *****************************



set.seed(123)
x<- matrix( runif( 30), 15,2)
Z<- matrix( rnorm(30), 15,2)
y<- rnorm( 15)*.01 + 5*(x[,1]**3 +
                           (x[,2]-.5)**2) +  (Z[,1] +Z[,2])*.001
knots<- x[1:5,]
weights<- runif(15)*10
y[5] <- y[5] + 3 # avoids GCV warning 

# compare to 
Krig( x,y,Z=Z, knots=knots, cov.function=Exp.cov,weights=weights,
verbose=FALSE, give.warnings=FALSE)-> out.new

Krig( x,y,Z=Z, knots=knots, cov.function=Exp.cov,weights=weights, 
          lambda=1)-> out.new2

# compare to each other
Krig.coef( out.new, lambda=1)-> look
# test for d coefficients
test.for.zero( out.new2$d, look$d, tag=" knots/weights fixed/default d coef")
# test for c coefficents
test.for.zero( out.new2$c, look$c, tag="knots/weights fixed/default c coef" )


# compute test using linear algebra

K<- Exp.cov( knots, knots)

T<- cbind( rep(1,15), x, Z)
X<- cbind( T, Exp.cov( x, knots))
lambda<-1.0
NN<- ncol( X)
H<- matrix( 0, NN, NN)
H[(1:5)+5, (1:5)+5] <- K

c(   solve(t(X)%*%(weights*X) + lambda*H)%*% t(X)%*% (weights*y) )-> temp
temp.c<- temp[6:10]
temp.d<- temp[1:5]

# test for d coefficients
test.for.zero( out.new2$d, temp.d, tag=" knots d coef")
# test for c coefficents
test.for.zero( out.new2$c, temp.c, tag="knots c coef" )


####### testing predict function 
hold1<- predict( out.new2, x=x, Z=Z, y=y)
hold2<- predict( out.new2, x=x, Z=Z, just.fixed=TRUE,y=y)
hold3<- predict( out.new2, x=x, Z=Z, drop.Z=TRUE,y=y)
hold4<- predict( out.new2, x=x, Z=Z, drop.Z=TRUE, just.fixed=TRUE,y=y)


hold<- X%*% temp
#  X%*% temp -  X[,4:5]%*% temp[c(4,5)]

test.for.zero( hold, hold1, tag="knots predict for null" )

hold<-T%*%temp.d
test.for.zero( hold, hold2, tag="knots predict for null" )

hold<-X%*%temp - X[,4:5] %*% temp[4:5]
test.for.zero( hold, hold3, tag="knots predict w/o Z" )

hold<-T[,1:3]%*%temp.d[1:3]
test.for.zero( hold, hold4, tag="knots predict for drift" )

######tests where coefficients  are recomputed from object
hold1<- predict( out.new,y=y, lambda=1.0,  x=x, Z=Z)
hold2<- predict( out.new,y=y, lambda=1.0, x=x, Z=Z, just.fixed=TRUE)
hold3<- predict( out.new, y=y, lambda=1.0, x=x, Z=Z, drop.Z=TRUE)
hold4<- predict( out.new, y=y, lambda=1.0, x=x, Z=Z, 
                      drop.Z=TRUE, just.fixed=TRUE)

hold<-X%*%temp
test.for.zero( hold, hold1, tag="predict for null" )

hold<-T%*%temp.d
test.for.zero( hold, hold2, tag="predict for null" )

hold<-X[,1:3] %*%temp.d[1:3] + X[,6:10] %*% temp.c
test.for.zero( hold, hold3, tag="predict for null" )

hold<-T[,1:3]%*%temp.d[1:3]               
test.for.zero( hold, hold4, tag="predict for null" )

options( echo=TRUE)
