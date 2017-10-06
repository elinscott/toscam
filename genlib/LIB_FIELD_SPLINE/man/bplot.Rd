\name{bplot}
\alias{bplot}
\title{
  boxplot  
}
\description{
Plots boxplots of several groups of data
and allows for placement at different horizontal or vertical positions or 
colors. It
is also flexible in the input object, accepting either a list or matrix.
}
\usage{
bplot(x, by,style = "tukey", outlier = TRUE, plot = TRUE, ...)
}
\arguments{
 \item{x}{
Vector, matrix,  list or data frame. A vector may be divided according
to the by argument. Matrices and data frames are separated by columns and
lists by components.   
}
\item{by}{
If x is a vector, an optional vector (either character or numerical)
specifying the categories to divide x into separate data sets.
}
\item{style}{
Type of boxplot default is "tukey". The other choice is "quantile" where
the whiskers are drawn to the 5 and 95 percentiles instead being based on
the inner fences. 
 
}
\item{outlier}{
If true outliers (points beyond outer fences) will be added to the plots. 
}
\item{plot}{
If false just returns a list with the statistics used for plotting the
box plots. 
}
\item{\dots}{
Other arguments controlling the boxplots ( passed to bplot.obj) these are
listed below. Other graphical plotting arguments not matched (e.g. yaxt)
are used in the call to plot to set up the initial plot if add=T. 
 
 \describe{
  \item{pos}{
  The boxplots will be plotted vertically and pos gives the x or y locations
  for their centers. If omitted the boxes are equally spaced at integer  
  values.  
  }
  
   \item{width}{
   Width of boxplots (in user coordinates) if omitted then the width is a
   reasonable fraction of the distance between boxes and is set by the 
   space argument. }
  
  \item{labels}{
  Labels under each boxplot. If missing the columns names or components of  
  x are used.}
  \item{las}{
  Orient the axis labels for bplot groups. Default is to 
  put them horizontal (las=1) if number of groups is less than 7 
  otherwise make them 
  perpendicular (las=2) to keep the labels from running into each other. 
  See help(par) for more about this option.}

  \item{add}{
   If true, do not create a new plots just add the boxplots to a current
   plot. Note that the pos argument may be useful in this case and should 
   be in the user coordinates of the parent plot.}

  \item{space}{Space between boxplots.}

  \item{sort.names}{If true plot the boxplot data set names are sorted in 
   alphabetic order by their labels.}

  \item{xlab}{Label for the x-axis}

  \item{ylab}{ Label for the y-axis}

  \item{label.cex}{
  Boxplot label  size where 1.0 is  normal size characters.
  If zero labels will not be added. }

  \item{xaxt}{
  Plotting parameter for x-axis generation. Default is not to produce
  an x-axis.}

  \item{horizontal}{
    If true draw boxplots horizontally the default is false, produce
    vertical box plots. }

  \item{lwd}{  Width(s) of lines in box plots. }

  \item{col}{ Color(s) of bplots.} 
  }
 }
}
\details{
This function was created as a complement to the usual S function for  
boxplots. The current function makes it possible to put the boxplots at  
unequal x or y positions. This is useful for visually grouping a large
set of  
boxplots into several groups. Also placement of the boxplots with  
respect to the axis can add information to the plot.  Another aspect  
is the emphasis on data structures for groups of data. One useful  
feature is the by option to break up the x vector into distinct  groups.   
If 5 or less observations are in a group the points themselves are  
plotted instead of a box.   

The function is broken into two steps: a call to stats.bplot to find the
statistics and a call to bplot.obj to plot the resulting object. 
The user is referred to describe.bplot to modify the statistics used and 
to draw.bplot.obj to modify how the bplot is drawn.

Finally to bin data into groups based on a continuous variable and to make
bplots of each group see bplot.xy. 
}
\seealso{
describe.bplot, draw.bplot.obj, stats.bplot, bplot.xy, bplot.obj
}
\examples{
#
set.seed(123)
temp<- matrix( rnorm(12*8), ncol=12)
pos<- c(1:6,9:14)
bplot(temp)
#
bplot( temp, pos=pos, labels=paste( "D",1:12), horizontal=TRUE)
#
# boxplots in red
bplot( temp, pos=pos, label.cex=0, horizontal=TRUE, col="red")
# add an axis
axis( 2)
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
