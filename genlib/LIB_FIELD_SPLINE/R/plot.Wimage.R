"plot.Wimage" <-
function (x, cut.min, graphics.reset = TRUE, common.range = FALSE, 
    color.table = tim.colors(128), Nlevel = NULL, with.lines = FALSE, 
    omd.width = 0.2, ...) 
{
    m <- nrow(x)
    n <- ncol(x)
    info <- Wimage.info(m, n, cut.min)
    if (is.null(Nlevel)) {
        Nlevel <- info$Lmax
    }
    old.par <- par(no.readonly = TRUE)
    par(mar = c(0, 1, 1, 0), ...)
    omd.bottom <- par()$cxy[2]
    figs <- rbind(c(0, 1 - omd.width, omd.bottom, 1), c(1 - omd.width, 
        1, omd.bottom, 1))
    split.screen(figs)
    ind <- split.screen(c(Nlevel + 1, 3), screen = 1)
    zr.common <- range(x, na.rm = TRUE)
    S1 <- info$S[1]:info$S[2]
    S2 <- info$S[3]:info$S[4]
    M <- 1:info$L[1, 1]
    N <- 1:info$L[1, 2]
#
# add box function
    add.boxes <- function() {
        if (with.lines) {
            xline(c(0.5, M + 0.5), col = "white", lwd = 0.5)
            yline(c(0.5, N + 0.5), col = "white", lwd = 0.5)
        }
    }

    screen(ind[1])
    image(M, N, x[S1, S2], zlim = zr.common, xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "", col = color.table)
    add.boxes()
    box()

    screen(ind[3])
    image(1:m, 1:n, x, xaxt = "n", yaxt = "n", col = color.table, 
        xlab = "", ylab = "")
    box()
#
#  add in the multiresolution partitioning of 
#  image matrix. 

    if( with.lines){
     temp<- info$L
      temp<- rbind( temp, c(m,n))
     for( k in 1:info$Lmax){
        xt<- c(temp[k,1]+.5, temp[k,1]+.5)
        yt<- c(0,temp[k+1,2]+.5) 
        segments( xt[1], yt[1],xt[2],yt[2], col="white", lwd=.5)
        yt<- c(temp[k,2]+.5, temp[k,2]+.5)
        xt<- c(0,temp[k+1,1]+.5) 
        segments( xt[1], yt[1],xt[2],yt[2], col="white", lwd=.5)
     }
    }

    save.zr <- matrix(NA, Nlevel, 2)
    screen.off <- ind[3]
    KK <- 1
    for (lev in (1:Nlevel)) {
        M <- 1:info$L[lev, 1]
        N <- 1:info$L[lev, 2]
        H1 <- info$H[lev, 1]:info$H[lev, 2]
        H2 <- info$H[lev, 3]:info$H[lev, 4]
        V1 <- info$V[lev, 1]:info$V[lev, 2]
        V2 <- info$V[lev, 3]:info$V[lev, 4]
        D1 <- info$Di[lev, 1]:info$Di[lev, 2]
        D2 <- info$Di[lev, 3]:info$Di[lev, 4]

        if (common.range) {
            zr <- zr.common
        }
        else {
            zr <- range(c(x[H1, H2], x[V1, V2], x[D1, D2]), na.rm = TRUE)
            save.zr[lev, ] <- zr
        }
        screen(screen.off + KK)
        image(M, N, x[H1, H2], zlim = zr, xaxt = "n", yaxt = "n", 
            col = color.table, xlab = "", ylab = "")
        add.boxes()
        box()
        KK <- KK + 1
        screen(screen.off + KK)
        image(M, N, x[V1, V2], zlim = zr, xaxt = "n", yaxt = "n", 
            col = color.table, xlab = "", ylab = "")
        add.boxes()
        box()
        KK <- KK + 1
        screen(screen.off + KK)
        image(M, N, x[D1, D2], zlim = zr, xaxt = "n", yaxt = "n", 
            col = color.table, xlab = "", ylab = "")
        add.boxes()
        box()
        KK <- KK + 1
    }
    screen(2)

    if (common.range) {
        image.plot(zlim = zr, legend.only = TRUE, smallplot = c(0.2, 
            0.3, 0.1, 0.5), col = color.table, graphics.reset = TRUE)
    }
    else {
        ind2 <- split.screen(c(Nlevel + 1, 1), screen = 2)
        screen(ind2[1])
        image.plot(zlim = zr, legend.only = TRUE, smallplot = c(0.2, 
            0.3, 0.1, 0.9), col = color.table)
        for (lev in 1:Nlevel) {
            screen(ind2[lev + 1])
            image.plot(zlim = save.zr[lev, ], legend.only = TRUE, 
                col = color.table, smallplot = c(0.2, 0.3, 0.1, 
                  0.9), graphics.reset = TRUE)
        }
    }
    if (graphics.reset) {
        close.screen(all = TRUE)
        par(cex.axis = 1)
        par(old.par)
        invisible()
    }
    else {
        cat("Still in split.screen mode", fill = TRUE)
        cat("  Use: screen.close( all=TRUE) to close this panel \n      after adding to plots", 
            fill = TRUE)
        return(matrix(ind, ncol = 3, byrow = TRUE))
    }
}

