"plot.Wimage" <-
function (x, cut.min , graphics.reset = TRUE, common.range = FALSE,
color.table=tim.colors(128),Nlevel=NULL,with.lines=FALSE,omd.width=.2,...) 
{

m<- nrow( x)
n<- ncol(x)
info<-  Wimage.info(m,n, cut.min)
if( is.null(Nlevel)){
Nlevel<- info$Lmax}

    old.par <- par(no.readonly=TRUE)
    par(mar = c(0, 1, 1, 0),...)

# split the screen vertically into the main panel space and 
# a strip to put the legends. 
omd.bottom<- par()$cxy[2]
figs<- rbind( 
             c(0, 1-omd.width, omd.bottom, 1), 
             c(1-omd.width, 1, omd.bottom, 1))

split.screen(figs)

# now split first screen into the panel
split.screen(c(Nlevel+1,3), screen=1)-> ind

# draw smooth basis coefficients. 

zr.common <- range(x, na.rm=TRUE)

S1<- info$S[1]:info$S[2]
S2<- info$S[3]:info$S[4]
M<- 1:info$L[1,1]
N<- 1:info$L[1,2]

  add.boxes<- function(){
                         if( with.lines){
                          xline( c( .5, M+.5), col="white", lwd=.5)
                          yline( c( .5, N+.5), col="white",lwd=.5)
                         }
                        }

# move to first plot
# smooth basis functions
screen(ind[1])
image(M,N,x[S1,S2], zlim = zr.common, xaxt = "n", yaxt = "n", 
       xlab="", ylab="",col=color.table)
add.boxes()

# full image
screen(ind[3])
 image( 1:m,1:n,x,xaxt = "n", yaxt = "n",col=color.table,
  xlab="", ylab="")

# save the ranges for each level to draw the legend strips later 
save.zr<- matrix( NA, Nlevel, 2)

screen.off<- ind[3]
KK<-1 
for( lev in (1:Nlevel)){
   
   M<- 1:info$L[lev,1]
   N<- 1:info$L[lev,2]
   H1<- info$H[lev, 1]:info$H[lev, 2]
   H2<- info$H[lev, 3]:info$H[lev, 4]
   V1<- info$V[lev, 1]:info$V[lev, 2]
   V2<- info$V[lev, 3]:info$V[lev, 4]
   D1<- info$Di[lev, 1]:info$Di[lev, 2]
   D2<- info$Di[lev, 3]:info$Di[lev, 4]

 if( common.range){zr<- zr.common}
  else{ 
       zr <- range(
        c(x[H1, H2], x[V1,V2], x[D1,D2]), na.rm=TRUE)
      save.zr[lev,] <- zr
  }
   screen( screen.off+KK)
   image(M,N,x[H1, H2], 
    zlim = zr, xaxt = "n", yaxt = "n",col=color.table,
    xlab="", ylab="")
   add.boxes()
   KK<- KK+1
   
   screen( screen.off+KK)
   image(M,N,x[V1,V2], 
    zlim = zr, xaxt = "n", yaxt = "n",col=color.table,
    xlab="", ylab="")
   add.boxes()
   KK<- KK+1

   screen( screen.off+KK)
   image(M,N,x[D1, D2], 
    zlim = zr, xaxt = "n",yaxt = "n", col=color.table,
    xlab="", ylab="")
   add.boxes()
   KK<- KK+1
   }

screen(2)
if( common.range){
image.plot(zlim=zr, legend.only=TRUE, smallplot=c(.2,.3,.1,.5),
 col=color.table, graphics.reset=TRUE)
}
else{
split.screen( c( Nlevel+1, 1), screen=2)-> ind2
screen( ind2[1])
 image.plot(zlim=zr, legend.only=TRUE, smallplot=c(.2,.3,.1,.9), 
           col=color.table)

 for ( lev in 1:Nlevel){
screen( ind2[lev+1])
   image.plot(zlim=save.zr[lev,], legend.only=TRUE,col=color.table,
     smallplot=c(.2,.3,.1,.9), graphics.reset=TRUE)}
 }

# reset margins 


 if (graphics.reset) {
   close.screen( all=TRUE)
   par( cex.axis=1.0)
   par(old.par)
   invisible()}
 else{
   return( matrix( ind, ncol=3, byrow=TRUE))
 }
}
"W.i2s" <-
function(ind, m, cut.min){

W.info( m, cut.min=cut.min)-> info

# figure out to which level the index addresses. 
M<- rep( 0, length( ind))
for(  kk in 1:info$Lmax){
M<- ifelse( ind > info$offset[kk], kk,M)
}

 
# compute off for levels taking into account 
# father wavelets and indices that are too large

off<- rep( 0, length(ind))
M.good<- (M>0)& (M<=info$Lmax)
off[ M.good]<- info$offset[M[M.good]]
off[ind>m] <- NA

list( 
level=M,i=ind-off,
m=m, cut.min=cut.min)

}
"Wimage.i2s" <-
function(ind, m,n,cut.min) 
{
if( is.matrix( ind)){
tempI<- ind[,1]
tempJ<- ind[,2]}
else{
tempJ<- ceiling(ind/m) 
tempI<- ind- (ceiling(ind/m) -1)*m
}

#cat( tempI, tempJ, fill=T)

Wimage.info( m=m, n=n, cut.min=cut.min)-> info
W.i2s( tempI, m=m, cut.min=info$S[2])-> tempM
W.i2s( tempJ, m=n, cut.min=info$S[4])-> tempN


# find H V or Di. 
flavor<- rep( 1, length( ind))
flavor<- ifelse( tempM$level> tempN$level, 1,2)
flavor[ tempM$level== tempN$level] <- 3
level<- ifelse( tempM$level> tempN$level, tempM$level, tempN$level)
Sind<- level==0

# this change makes the array indexing work below
flavor[ Sind] <- NA
level[Sind] <- NA

I<- tempI - info$offset.m[cbind(level,flavor)]
J<- tempJ - info$offset.n[cbind(level,flavor)] 
#
#set bad values to NAs
I<- ifelse( (I<1)|(I>info$L[level,1]),NA, I)

#cat( level, J,info$L[level,2], fill=TRUE)
J<- ifelse( (J<1)|(J>info$L[level,2]),NA, J)


# reset level and flavor for fathers
flavor[ Sind] <- 0
level[ Sind] <- 0
I[Sind] <- tempI[Sind]
J[Sind] <- tempJ[Sind]

 
list(  i=I, j=J,
level=level, flavor=flavor,cut.min=cut.min, m=m,n=n)
}
"Wimage.info" <-
function (m=128,n=m, cut.min = 4)
{
    NN <- n
    MM <- m
    level <- 1
    while (min(c(NN, MM)) > cut.min) {
        level <- level + 1
        NN <- NN/2
        MM <- MM/2
    }
    nlevel<-level-1
#   print( nlevel)

    V<-Di<-H<- matrix(NA, nlevel, ncol=4)
    N<- matrix( NA, nlevel, ncol=2 )
    S<- matrix( NA, 1, ncol=4)
    offset.m<- matrix( NA, nlevel, ncol=3)
    offset.n<- matrix( NA, nlevel, ncol=3)
                                                                               
    n2 <- NN
    m2 <- MM
    n3 <- NN * 2
    m3 <- MM * 2
    n1 <- 1
    m1 <- 1                                       
    level <- 1
    
    S<- c( 1, m2, 1,n2)
    
    while (n3 <= n & m3 <= m) {
#     cat( level, fill=TRUE)
#     cat( m1,m2,m3, fill=TRUE)
#     cat( n1,n2,n3, fill=TRUE)
# H and V are defined so that an ordinary image plot 
# gives the expected orientationd for basis functions
#
        H[level,]<-c((m2 + 1),m3, n1,n2)
        V[level,]<-c(m1,m2, (n2 + 1),n3)
        Di[level,]<- c((m2+ 1),m3, (n2 + 1),n3)
        offset.m[level,] <- c(H[level,1], V[level,1], Di[level,1]) -1
        offset.n[level,] <- c(H[level,3], V[level,3], Di[level,3]) -1
                                                                                          
        N[level,] <- c( (m3-m2), c(n3-n2))
                            
        level <- level + 1
        n2 <- n3
        n3 <- n2 * 2
        m2 <- m3
        m3 <- m2 * 2
    }
list( m=m,n=n,cut.min=cut.min,S=S,H=H, V=V, Di=Di, L=N, Lmax =nlevel,
offset.m=offset.m,
offset.n=offset.n)
}
"Wimage.s2i" <-
function(i,j,level,flavor,m,n, cut.min, mat=TRUE){

Wimage.info(m=m,n=n,cut.min=cut.min)-> ind

if( is.list(i)){
j<- i$j
level<- i$level
flavor<- i$flavor
i<- i$i}


if( length(i)!= length(j) ) {
   stop( "i and j must be same length")}

if( length(level)==1){ 
   level <- rep( level, length(i))}
if( length(flavor)==1){ 
   flavor <- rep( flavor, length(i))}

#
# reset structure elements to NA that are out of bounds w/r to level
Sind<- level==0
level[(level > ind$Lmax)| Sind] <- NA
flavor[Sind] <- NA

hold<- i + ind$offset.m[cbind(level,flavor)] + 
           (j + ind$offset.n[cbind(level,flavor)]-1)*ind$m

# handle the father wavelet cases
hold[Sind] <-  i[Sind] +(j[Sind]-1)*ind$m
 

# reset structure elements to NA that are out of bounds w/r to position
hold[ (i<1)|(j<1)] <- NA
hold[ i> ind$L[level,1] ]<- NA
hold[j> ind$L[level,2] ]<- NA

hold

if( mat){ 
return(  cbind( hold- (ceiling(hold/m)-1)*m,ceiling(hold/m))  )
}
else{return( hold)}

}
"W.info" <-
function (m=128, cut.min = 4) 
{
    MM <- m
    level <- 1
    while (MM > cut.min) {
        level <- level + 1
        MM <- MM/2
    }
    nlevel<-level-1
#   print( nlevel)
    
    N<- matrix( NA, nlevel, ncol=1 )
    H<- matrix( NA, nlevel, ncol=2 )
    S<- matrix( NA, nrow=1, ncol=2)
    offset<- matrix( NA, nlevel, ncol=1)
 
    m2 <- MM
    m3 <- MM * 2
    m1 <- 1

    level <- 1

    S<- c( 1, m2) 
      off<- m2
    while (m3 <= m) {

        H[level,]<-c((m2+1),m3)
        offset[level,1] <- off 
        N[level,1] <- c( (m3-m2)) 
        
        off<- off+ m3-m2 
        level <- level + 1
        m2 <- m3
        m3 <- m2 * 2
    }
list( m=m,cut.min=cut.min,S=S,H=H, L=N, Lmax =nlevel, 
offset=offset)

}
"W.s2i" <-
function(i,level,m, cut.min){

# find the offsets in the array for different levels 
W.info(m,cut.min=cut.min)-> ind


# some checks 
if( length(level)==1){ 
   level <- rep( level, length(i))}

if( length(level)!=length( i))
stop("number of levels different from number of positions")

# indicator for smoothed (father) basis functions
indS<- level==0 

# reset level for fathers to one 
level[ indS] <- 1
# reset structure elements to NA that are out of bounds w/r to level
level[(level>ind$Lmax)] <- NA

# add postion in level to the level offset 

hold<- ifelse( indS, i, i + ind$offset[level] )

# catch any postions that exceed the allowed limits
hold[ (i>ind$L[level])|(i<1)] <- NA

hold

}
