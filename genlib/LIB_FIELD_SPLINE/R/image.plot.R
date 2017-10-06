"image.plot" <-
function (..., add=FALSE, nlevel = 64, 
    horizontal = FALSE, 
    legend.shrink = 0.9, 
    legend.width = 1.2, legend.mar = ifelse( horizontal, 3.1, 5.1),
    legend.lab=NULL, 
    graphics.reset = FALSE, 
    bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
    col = tim.colors(nlevel), lab.breaks=NULL, axis.args=NULL, legend.args=NULL) 
{

# Thanks to S. Koehler and  S. Woodhead
# for comments on making this a better function
#

# save current graphics settings
    old.par <- par(no.readonly = TRUE)

#  figure out zlim from passed arguments

    info <- image.plot.info(...)

    if (add) {
        big.plot <- old.par$plt
    }
  
    if (legend.only) {
        graphics.reset <- TRUE
    }
  
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1,5.1)
    }
#
# figure out how to divide up the plotting real estate. 
# 
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
#
# bigplot are plotting region coordinates for image
# smallplot are plotting coordinates for legend 
 
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot

#
# draw the image in bigplot, just call the R base function 
#
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot) }

        image(..., add = add, col = col)

        big.par <- par(no.readonly = TRUE)
    }

##
## check dimensions of smallplot

    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }

# Following code draws the legend using the image function 
# and a one column image. 
    
# calculate locations for colors on legend strip
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz) / nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints

    iz <- matrix(iy, nrow = 1, ncol = length(iy))

# extract the breaks
    breaks<- list(...)$breaks
   
# draw either horizontal or vertical legends. 
# using either suggested breaks or not -- a total of four cases. 
# 
   
# next par call sets up a new plotting region just for the legend strip
# at the plot coordinates 
         par(new = TRUE, pty = "m", plt = smallplot, err = -1)

                                        
# create the argument list to draw the axis
#  this avoids 4 separate calls to axis and allows passing extra
# arguments.

    if( is.null( breaks)){
           axis.args<- c( list( side= ifelse(horizontal,1,4),
                            mgp = c(3, 1, 0),
                            las = ifelse(horizontal, 0,2)),
                           axis.args)}
    else{
# add axis but only put label where there are breaks
# default labels for breaks are just the numeric values
                              
           if( is.null(lab.breaks)){ lab.breaks<- format( breaks)}

           axis.args<- c( list( side= ifelse(horizontal,1,4),
                                  mgp = c(3, 1, 0),
                                 las = ifelse(horizontal, 0,2),
                                 at=breaks,
                                 labels=lab.breaks) ,
                          axis.args)}    
# 
# draw color scales the four cases are horizontal/vertical breaks/no breaks
# add a label if this is passed.
    
    if (!horizontal) {
    
          if( is.null( breaks)){
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)}
          else{
            image(ix, iy, iz, 
                 xaxt = "n", yaxt = "n", xlab = "",  ylab = "", 
                 col = col, breaks=breaks )}
         }
    else {
      
          if( is.null( breaks)){
             image(iy, ix, t(iz), xaxt="n",yaxt = "n", xlab = "", ylab = "", 
                 col = col)}
          else{
             image(iy, ix, t(iz), 
                 xaxt = "n", yaxt = "n", xlab = "",  ylab = "", 
                 col = col, breaks=breaks )}
         }
    

# 
#
# now add the axis to the legend strip.
             
        do.call( "axis", axis.args)


# add a box around legend strip
        box() 
#
#
# add a label to the axis if information has been  supplied
# using the mtext function. The arguments to mtext are
# passed as a list like the drill for axis (see above)
#
            if( !is.null( legend.lab) ){         
             legend.args<-list( 
                    text= legend.lab,side= ifelse(horizontal, 1, 4), 
                    line=legend.mar - 2)
#                    just guessing at a good default for line! 
            }
#
# add the label using mtext function             

            if( !is.null( legend.args)){
            do.call( mtext, legend.args) }
#
#
# clean up graphics device settings
# reset to larger plot region with right user coordinates. 

    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}

