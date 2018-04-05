"image.plot.plt" <-
function (x, add = FALSE, legend.shrink = 0.9, legend.width = 1, 
    horizontal = FALSE, legend.mar = NULL, bigplot = NULL, smallplot = NULL, 
    ...) 
{
    old.par <- par(no.readonly = TRUE)
    if (is.null(smallplot)) 
        stick <- TRUE
    else stick <- FALSE
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
# compute how big a text character is 
    char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2], 
        par()$cin[1]/par()$din[1])

# This is how much space to work with based on setting the margins in the 
# high level par command to leave between strip and big plot
    offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])

# this is the width of the legned strip itself.
    legend.width <- char.size * legend.width
# this is room for legend axis labels
    legend.mar <- legend.mar * char.size

# smallplot is the plotting region for the legend. 
    if (is.null(smallplot)) {
        smallplot <- old.par$plt
        if (horizontal) {
            smallplot[3] <- legend.mar
            smallplot[4] <- legend.width + smallplot[3]
            pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
            smallplot[1] <- smallplot[1] + pr
            smallplot[2] <- smallplot[2] - pr
        }
        else {
            smallplot[2] <- 1 - legend.mar
            smallplot[1] <- smallplot[2] - legend.width
       
            pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
            smallplot[4] <- smallplot[4] - pr
            smallplot[3] <- smallplot[3] + pr
        }
    }
    if (is.null(bigplot)) {
        bigplot <- old.par$plt
        if (!horizontal) {
            bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
        }
        else {
            bottom.space <- old.par$mar[1] * char.size
            bigplot[3] <- smallplot[4] + offset
        }
    }
    if (stick & (!horizontal)) {
        dp <- smallplot[2] - smallplot[1]
        smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
        smallplot[2] <- smallplot[1] + dp
    }
    return(list(smallplot = smallplot, bigplot = bigplot))
}

