"bplot.obj" <-
function (data, pos = NA, width = NULL, labels = NULL, las = NULL, 
    add = FALSE, space = 0.25, sort.names = FALSE, xlab = "", 
    ylab = "", label.cex = 1, xaxt = "n", outlier = TRUE, horizontal = FALSE, 
    lwd = NA, col = NA, ...) 
{
    cols <- length(data)
    if (is.na(lwd)) {
        lwd <- par()$lwd
    }
    if (is.na(col)) {
        col <- par()$col
    }
    if (length(col) == 1) {
        col <- rep(col, cols)
    }
    if (length(lwd) == 1) {
        lwd <- rep(lwd, cols)
    }

    range.data <- c(NA, NA)

    if (is.null(labels)) {
        labels <- names(data)
    }

    if (is.na(pos[1])) {
        pos <- 1:cols
        if (sort.names) {
            pos <- order(labels)
        }
    }

    if (is.null(width)) {
        width <- min(diff(sort(pos))) * space
        if (cols == 1) 
            width <- space
    }

    if (length(width) == 1) 
        width <- rep(width, cols)
# determine limits for plotting bplots
    if (!add) {
        for (k in 1:cols) {
# if plotting oultiers too use the full range
# outliers are not plotted then just use limits of boxplot whiskers (in $bb)
         if( outlier){
            range.data <- range(c(range.data, data[[k]]$range), na.rm = TRUE)}
         else{
            range.data <- range(c(range.data, data[[k]]$bb), na.rm = TRUE)}


        }
        temp1 <- range.data
        temp2 <- range(c(pos - (0.5 * width)/space, pos + (0.5 * 
            width)/space))
        if (horizontal) {
            plot(temp1, temp2, type = "n", yaxt = xaxt, xlab = xlab, 
                ylab = ylab, ...)
        }
        else {
            plot(temp2, temp1, type = "n", xaxt = xaxt, xlab = xlab, 
                ylab = ylab, ...)
        }
    }
    for (i in 1:cols) {
        draw.bplot.obj(data[[i]], width[i], pos[i], outlier = outlier, 
            horizontal = horizontal, lwd = lwd[i], col = col[i])
    }
    if (label.cex > 0) {
        if (is.null(las)) {
            if (length(labels) > 7) {
                las <- 2
            }
            else {
                las <- 1
            }
        }
        if (horizontal) {
            axis.loc <- 2
        }
        else {
            axis.loc <- 1
        }
        axis(axis.loc, pos, labels, tick = FALSE, las = las, 
            adj = 0.5, cex = label.cex)
    }
    invisible()
}

