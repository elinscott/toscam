"Wimage.info.plot" <-
function( m,n,cut.min){
temp<- Wimage.info( m,n, cut.min)
Nmax<- max( m,n)
plot( c(1,Nmax), c(1,Nmax), type="n")
rect( temp$S[1], temp$S[3], temp$S[2], temp$S[4], col=1)
for (  k in 1: temp$Lmax){

MM<- temp$L[k,1]
NN<- temp$L[k,2]

rect( temp$H[k,1], temp$H[k,3], temp$H[k,2], temp$H[k,4], col=3)
rect( temp$V[k,1], temp$V[k,3], temp$V[k,2], temp$V[k,4], col=4)
rect( temp$Di[k,1], temp$Di[k,3], temp$Di[k,2], temp$Di[k,4], col=2)

# add to H
m0<- temp$offset.m[k,1]
n0<- temp$offset.n[k,1]
points( m0+c(1,MM), n0+c(1,NN), pch="O",cex=1)

# add to V
m0<- temp$offset.m[k,2]
n0<- temp$offset.n[k,2]
points( m0+c(1,MM), n0+c(1,NN), pch="O",cex=1)

# add to Di
m0<- temp$offset.m[k,3]
n0<- temp$offset.n[k,3]
points( m0+c(1,MM), n0+c(1,NN), pch="O",cex=1)
}
}
