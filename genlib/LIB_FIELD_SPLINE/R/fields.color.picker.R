
fields.color.picker<-function(){

  c( mar=c(0,0,3,0) )

# names of colors in default graphics options. 
 clab<- colors()
 n<- length( clab)
 N<- ceiling( sqrt(n) )
 M<- N
 temp<- rep( NA,M*N)
 temp[1:n] <- 1:n
 z<- matrix(temp, M,N)

# matrix of all colors

 image(seq(.5,M+.5,,M+1), seq(.5,N+.5,,N+1)
       , z,  col=clab, axes=FALSE, xlab="", ylab="")

  cat("Use mouse to identify color",fill=TRUE)
  locator( 1)-> loc; i<- round( loc$x);j<- round( loc$y);
  ind<- z[i,j] 

 points( i,j,  col= clab[ind], cex=4, pch="O")
 points( i,j, pch="+", col= "black", cex=1)

 mtext( side=3, text=clab[ind], col=clab[ind], line=1, cex=2)

# write out RGB values to console

 cat("ID ", ind, " name ", clab[ind],fill=TRUE)
 cat( "RGB", col2rgb(clab[ind])/256, fill=TRUE)
 temp<- signif(col2rgb(clab[ind])/256, 3)

# This line is  marginally in  LaTeX format to define color
 cat( clab[ind],
  " {rgb}{", temp[1], 
         ",", temp[2], 
         ",", temp[3], 
          "}", fill=TRUE )

}


