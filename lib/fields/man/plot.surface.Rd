\name{plot.surface}
\alias{plot.surface}
\title{
  Plots a surface 
}
\description{
Plots a surface object in several different ways to give 3-d
information e.g. a contour plots, perspective plots.  
}
\usage{
plot.surface(x, main = NULL, type = "C", zlab = NULL, xlab = NULL,
    ylab = NULL, levels = NULL, zlim = NULL, graphics.reset = NULL,
     labcex = 0.6, add.legend=TRUE, ...)
}
\arguments{
\item{x}{
A surface object. At the minimum a list with components x,y and z 
in the same form as the input list for the standard contour, persp
or image functions. This can also be an object from predict.surface. 
}
\item{main}{
Title for plot. 
 
}
\item{type}{
type="p" for a perspective/drape plot, 
type="I" for an image plot with a legend
strip (see image.plot). type="C" is the "I" option with contours lines
added. type="b" gives both  "p" and "C" as a 2X1 panel }
\item{zlab}{
z-axes label 
}
\item{xlab}{
x-axes label 
}
\item{ylab}{
y-axes labels 
}
\item{levels}{
Vector of levels to be passed to contour function. 
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting. 
Default is to reset if type ="b" but not for the single plot options.  
}
\item{zlim}{
Sets z limits on perspective plot. }
\item{labcex}{
Label sizes for axis labeling etc. 
}
\item{add.legend}{ If TRUE adds a legend to the draped perspective plot}
\item{\dots}{
Other graphical parameters that are passed along to either drape.persp 
or image.plot 
}

}
\seealso{
surface, predict.surface, as.surface, drape.plot, image.plot  
}
\examples{
fit<- Tps( BD[,1:4], BD$lnya) # fit surface to data 

# surface of variables 2 and 3 
#    holding 1 and 4 fixed at their median levels
 out.p<-predict.surface(fit, xy=c(2,3))  
 plot.surface(out.p, type="C") # surface plot  

}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
