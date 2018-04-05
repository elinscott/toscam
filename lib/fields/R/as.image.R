"as.image" <-
function (Z, ind = NULL, grid = NULL, x = NULL, nrow = 64, ncol = 64, 
    weights = NULL, na.rm = FALSE,nx=NULL,ny=NULL, boundary.grid=FALSE) 
{

# NOTE that throughout ind is a two column integer matrix of 
# discretized locations in the image matrix. 
# from different reductions due to NAs the final 
# versions of Z, x, and weights may be a subset of the 
# passed versions. 

# Thanks to J. Rougier for fixing bugs in this function.   
  
# set some default values for arguments
#
    if (!is.null(ind)) 
        x <- ind

    if( is.null( weights)) {
      weights<- rep( 1, length(Z))}
#
# use values of nx ny if passed; these are just different names for 
#  nrow, ncol.
 
    if( !is.null(nx)) nrow<- nx
    if( !is.null(ny)) ncol<- ny

##### end of setting defaults

    
# check for missing values in Z and  na.rm==FALSE 
    if( any( is.na(Z) ) & !na.rm ) {
           stop("missing values in Z, set na.rm=TRUE")}

#
#  if there are missing overwrite Z, x, ind and weights. 
#
      if (any( is.na(Z)) ) {
         Z.good <- !is.na(Z)
         Z <- Z[Z.good]
         x <- x[Z.good, ]
         ind<- ind[Z.good,]
         weights <- weights[Z.good]
      }
#
# check for x or weights having missing values 
# we do not like these ...

     if( any(is.na(weights)) | any( is.na( c( x)))) { 
           stop("missing values in weights or x")}


# discretize locations to grid boxes
# this function will also create a default grid based on range of
# locations if is NULL
#
        temp <- discretize.image(x,m=nrow, n=ncol,
                                grid = grid, boundary.grid=boundary.grid)
    
        nrow<- temp$m
        ncol<- temp$n
        ind <- temp$index
        grid<- temp$grid
    
# ind is a two column matrix with the index of the x and y grid points. 
# NOTE points outside of grid are NAs.
#
#

# if any of the ind's rows are NA's it means that the x's were 
# outside the range of the grid. 
# pare down arguments and eliminate these points. 
# 
   good.ind<- !( is.na(ind[,1]) | is.na( ind[,2]) )
   if( any( !good.ind) ){
           warning( "Some locations are outside the grid limits") 
           ind<- ind[good.ind,]
           Z<- Z[good.ind]
           weights<- weights[good.ind] }

# 
# find unique set of boxes for the locations 
#
# 
      rep.info <- cat.matrix(ind)
    uniquerows <- !duplicated(rep.info)

#
# compute weighted means where there are replicates 
#
# NOTE that in the assignments below we use the
# fact that a 2 column matrix (i.e. ind) is interpreted as a multiple index. 
#

    if (sum(uniquerows) < length(Z)) {
# this means that some Z's are in the same box 
        ind <- ind[uniquerows, ]
        temp <- fast.1way(rep.info, Z, w = weights)

# over write Z with  weighted means
        Z <- temp$means
        Ncell <- temp$n
        temp2 <- matrix(0, nrow = nrow, ncol = ncol)
        temp2[ind] <- Ncell
        temp3 <- matrix(NA, nrow = nrow, ncol = ncol)
        temp3[ind] <- temp$w.means
    }
    else {
        temp2 <- matrix(0, nrow = nrow, ncol = ncol)
        temp2[ind] <- 1
        temp3 <- matrix(NA, nrow = nrow, ncol = ncol)
        temp3[ind] <- 1
    }


    temp <- matrix(NA, nrow = nrow, ncol = ncol)
    temp[ind] <- Z

    call <- match.call()
    list(x = grid$x, y = grid$y, z = temp, call = call, ind = ind, 
        N = temp2, weights = temp3)
}

