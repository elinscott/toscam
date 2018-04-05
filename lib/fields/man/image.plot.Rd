\name{image.plot}
\alias{image.plot}
\title{
  Draws image plot with a legend strip for the color scale.  
}
\description{
This function combines the R image function with some automatic 
placement of a legend. This is done by splitting the plotting region  
into two parts. Putting the image in one and the legend in the other.  
}
\usage{

image.plot(..., 
                 add = FALSE, nlevel = 64, horizontal = FALSE,
                 legend.shrink = 0.9, legend.width = 1.2, legend.mar =
                 ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
                 graphics.reset = FALSE, bigplot = NULL, smallplot =
                 NULL, legend.only = FALSE, col = tim.colors(nlevel),
                 lab.breaks = NULL, axis.args = NULL, legend.args =
                 NULL)}

\arguments{ \item{\dots}{
The usual arguments to the image function. This includes the use of 
the breaks argument for an unequal color scale. 
}
\item{add}{
If true add image and a legend strip to the existing plot. 
}
\item{nlevel}{
Number of color levels used in legend strip  
}
\item{legend.shrink}{
 Amount to shrink the size of legend relative to the full height or width
of the plot.  
}
\item{legend.width}{
Width in characters of the legend strip. Default is 1.2, a little bigger 
that the width of a character. }

\item{legend.mar}{
Width in characters of legend margin that has the axis. Default is 5.1 
for a vertical legend and 3.1 for a horizontal legend.}
\item{legend.lab}{ Label for the axis of the color legend. Default is no
label as this is usual evident from the plot title.}

\item{graphics.reset}{
 If
FALSE (default)  the plotting region ( plt in par) will not be reset 
and
one can add more information onto the image plot. (e.g. using functions
such as points or lines.)  If TRUE will reset plot parameters to the
values before entering the function.  
}
\item{horizontal}{
If false (default) legend will be a vertical strip on the right side. If
true the legend strip will be along the bottom.  
}
\item{bigplot}{
Plot coordinates for image plot. If not passed
these will
be determined within the function.  
}
\item{smallplot }{
Plot coordinates for legend. If not passed these will be determined within
the function.  
}
\item{legend.only}{
If TRUE just add the
legend to a the plot in the plot region defined by the coordinates in
smallplot. In the absence of other information the range for the legend 
is determined from the \code{zlim} argument. 
}
\item{col}{
Color table to use for image ( see help file on image for details).
Default is a pleasing range of 64 divisions suggested by Tim Hoar and is similar to 
the MATLAB (TM) jet color scheme.
}

\item{lab.breaks}{ If breaks are supplied these are text string labels
to put at each break value. This is intended to label axis on a
transformed scale such as logs.}
\item{axis.args}{Additional arguments for the axis function used to 
create the legend axis. (See example below adding a log scaling.)}

\item{legend.args}{Arguments for a complete specification of the 
legend label. This is in the form of list and is just passed to the 
mtext function. Usually this will not be needed. 
(See example below.)}

}
\section{Side Effects}{
After exiting, the
plotting region may be changed to make it possible to add more features to
the plot. To be explicit, \code{par()\$plt} may be changed to reflect a 
smaller plotting region that has accommodated room for the legend subplot.  
}
\details{
 If the z component is a matrix then the user should be aware that this
function locates the matrix element z[i,j] at the grid locations (x[i],
y[j]) this is very different than simply listing out the matrix in the
usual row column tabular form. See the example below for more details of
this difference in formatting. What does one do if you don't really have
the "z" values on a regular grid? See the functions \code{quilt.plot.Rd}
and \code{as.image} to discretize irregular observations to a grid. 


It is surprising how hard it is just to automatically add the
legend! All "plotting coordinates" mentioned here are in device
coordinates. The plot region is assumed to be [0,1]X[0,1] and plotting
regions are defined as rectangles within this square. We found these
easier to work with than user coordinates.  

\code{legend.width} and \code{legend.mar} are in character spaces. These
units are helpful in thinking about axis labels that will be put into
these areas. To add more or less space between the legend and the image
plot alter the mar parameters. The default mar settings (5.1,5.1,5.1,2.1)
leaves 2.1 spaces for vertical legends and 5.1 spaces for horizontal
legends. Changing the plot margins directly replaces the offset 
argument in the older version of this function. 

There are always problems with
default solutions to placing information on graphs but the choices made
here may be useful for most cases. The most annoying thing is that after
using plot.image and adding information the next plot that is made may
have the slightly smaller plotting region set by the image plotting.
The user should set \code{reset.graphics=TRUE} to avoid the plotting size 
some from changing. The disadvantage, however, of resetting the graphics 
is that one can no longer add additional graphics elements to the image 
plot. Note that filled.contour always resets the graphics but provides 
another mechanism to pass through plotting commands. Apparently 
filled.contour, while very pretty, does not work for multiple plots. 
\code{levelplot} that is part of the lattice package has a very 
similar function to image.plot and a formula syntax in the call. 

The strategy for \code{image.plot} is simple, divide the plotting region
into two smaller regions \code{bigplot} and \code{smallplot}. The image
goes in one and the legend in the other. This way there is always room for
the legend. Some adjustments are made to this rule by not shrinking the
\code{bigplot} if there is already room for the legend strip and also
sticking the legend strip close to the image plot. One can specify the
plot regions explicitly by \code{bigplot} and \ \code{smallplot} if the
default choices do not work. There may be problems with small plotting
regions in fitting both of these elements in the plot region and one may 
have to change the default character sizes or margins to make things fit.

By keeping the zlim argument the same across images one can generate the
same color scale. (See the image help file.)  One useful technique for a
panel of images is to just draw the images with \code{image} 
and then use image.plot to add a legend to the last plot. (See example 
below for messing with the outer margins to make this work.) 
Usually a square plot (pty="s") done in a rectangular plot region will
have room for the legend stuck to the right side without any other 
adjustments. See the example below and the code for \code{plot.Wimage} 
for more complicated arrangements of multiple image plots and summary 
legends.  

Note that to add just the legend strip all the numerical information one 
needs is the zlim argument! 
We like \code{tim.colors} as a default color scale. 
The the topographic color scale 
(\code{topo.colors}) is also a close second showing our geophysical basis. 
See also \code{terrain.colors} for a subset. For using other color choices 
see how the \code{nlevels} argument figures into the legend and main plot
number of colors. 

}

\seealso{ image,filled.contour, quilt.plot, plot.surface, 
add.image, colorbar.plot, tim.colors } 

\examples{
x<- 1:10; y<- 1:15; z<- outer( x,y,"+") 


image.plot(x,y,z) 
# or obj<- list( x=x,y=y,z=z); image.plot(obj)
# now add some points on diagonal with some clipping anticipated 
   points( 5:12, 5:12, pch="X", cex=3)

image.plot(x,y,z, legend.lab="inches")

# adding breaks and distinct colors
#        (lab.breaks are optional)

brk<- quantile( c(z))
image.plot(x,y,z, breaks=brk, col=rainbow(4), 
             lab.breaks=names(brk))
#
# compare to 

quantile(c(z), c( .05, .1,.5, .9,.95))-> zp

image.plot(x,y,z, 
   axis.args=list( at=zp, labels=names(zp) ) )


# a log scaling for the colors
ticks<- c( 1, 2,4,8,16,32)
image.plot(x,y,log(z), axis.args=list( at=log(ticks), labels=ticks))


#
#fat (5 characters wide) and short (50\% of figure)  color bar on the bottom
   image.plot( x,y,z,legend.width=5, legend.shrink=.5, horizontal=TRUE) 

# adding label with all kinds of additional arguments.
# use side=4 for vertical legend and side= 1 for horizontal legend
# to be paralell to axes. See help(mtext).

image.plot(x,y,z, 
       legend.args=list( text="unknown units",
     col="magenta", cex=1.5, side=4, line=2))



set.panel()

# Here is quick but quirky way to add a common legend to several plots. 
# The idea is leave some room in the margin and then over plot in this margin

par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 2,2) # 2X2 matrix of plots

# now draw all your plots using usual image command
for (  k in 1:4){
image( matrix( rnorm(150), 10,15), zlim=c(-4,4), col=tim.colors())
}

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=c(-4,4)) 

# image.plot tricked into  plotting in margin of old setting 

set.panel() # reset plotting device

#
# Here is a more learned strategy to add a common legend to a panel of
# plots  consult the split.screen help file for more explanations.
# For this example we  draw two
# images top and bottom and add a single legend color bar on the right side 



# first divide screen into the figure region and legend colorbar on the 
# right to put a legend. 

   split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))

# now divide up the figure region 
   split.screen(c(2,1), screen=1)-> ind

zr<- range( 2,35)
# first image
   screen( ind[1])
   image( x,y,z, col=tim.colors(), zlim=zr)

# second image
   screen( ind[2])
   image( x,y,z+10, col=tim.colors(), zlim =zr)

# move to skinny region on right and draw the legend strip 
   screen( 2)
   image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
   col=tim.colors())

   close.screen( all=TRUE)


# you can always add a legend arbitrarily to any plot;
# note that here the plot is too big for the vertical strip but the
# horizontal fits nicely.
plot( 1:10, 1:10)
image.plot( zlim=c(0,25), legend.only=TRUE)
image.plot( zlim=c(0,25), legend.only=TRUE, horizontal =TRUE)

# combining the  usual image function and adding a legend
# first change margin for some more room
\dontrun{
par( mar=c(10,5,5,5))
image( x,y,z, col=topo.colors(64))
image.plot( zlim=c(0,25), nlevel=64,legend.only=TRUE, horizontal=TRUE,
col=topo.colors(64))
}
#
# 
# sorting out the difference in formatting between matrix storage 
# and the image plot depiction

A<- matrix( 1:48, ncol=6)
# Note that matrix(c(A), ncol=6) == A
image.plot(1:8, 1:6, A)
# add labels to each box 
text( c( row(A)), c( col(A)), A)
# and the indices ...
text( c( row(A)), c( col(A))-.25,  
   paste( "(", c(row(A)), ",",c(col(A)),")", sep=""), col="grey")

# "columns" of A are horizontal and rows are ordered from bottom to top!
#
# matrix in its usual tabular form where the rows are y  and columns are x

image.plot( t( A[6:1,]), axes=FALSE)



 
}
\keyword{hplot}
% docclass is function
