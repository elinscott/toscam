"Krig.replicates" <-
function(out, verbose = FALSE)
{
        rep.info <- cat.matrix(out$x)

        if(verbose){
                cat( "replication info", fill=TRUE)
                print(rep.info)}

        uniquerows <- !duplicated(rep.info)
        if(sum(uniquerows) == out$N) {
                shat.rep <- NA
                shat.pure.error <- NA
                pure.ss <- 0
                yM <- out$y
                weightsM <- out$weights
                xM <- as.matrix(out$x[uniquerows,  ])
                ZM<- out$Z
        }
        else {
                rep.info.aov <- fast.1way(rep.info, out$y, out$weights)
                shat.pure.error <- sqrt(rep.info.aov$MSE)
                shat.rep <- shat.pure.error
                yM <- rep.info.aov$means
                weightsM <- rep.info.aov$w.means
                xM <- as.matrix(out$x[uniquerows,  ])

# choose some Z's for replicate group means
                if( !is.null(out$Z)){
                ZM<-  as.matrix(out$Z[ uniquerows,]) }
                else{ ZM<- NULL}

                pure.ss <- rep.info.aov$SSE
                if(verbose)
                        print(rep.info.aov)
        }
        return(list(yM=yM, xM=xM,ZM= ZM,
                weightsM=weightsM, uniquerows=uniquerows,
                shat.rep=shat.rep, shat.pure.error=shat.pure.error, pure.ss=pure.ss,
                 rep.info=rep.info))
}

