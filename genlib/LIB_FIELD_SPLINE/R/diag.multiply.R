
# create a special matrix multiply for diagonal
# matrices. Diagonal matrix assumed to be just a vector.
# see tests directory for tests
# NOTE: this is not a symmetric operation:
#  when a left vector is given it is a diagonal matrix
#  when a right vector is given it is a vector. 
#

setGeneric("%d*%",function(x,y,...)standardGeneric("%d*%"))

setMethod("%d*%",signature(x="matrix",y="matrix"),
function(x,y){x%*%y} )

setMethod("%d*%",signature(x="matrix",y="numeric"),
function(x,y){ x%*%y} )

setMethod("%d*%",signature(x="numeric",y="matrix"),
function(x,y){x*y} )

setMethod("%d*%",signature(x="numeric",y="numeric"),
function(x,y){cbind(x*y)} )
