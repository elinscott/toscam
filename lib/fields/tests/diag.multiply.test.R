library( fields)
options( echo=FALSE)
set.seed( 234)
n <- 5
m <- 4
mat <- array(rnorm(n*m),c(n,m))
mat2 <- array(rnorm(n*m),c(m,n))
vec <- rnorm(n)
vec2 <- rnorm(n)

test.for.zero( mat2 %*% mat, mat2%d*%mat, tol=1e-8 )

test.for.zero( (diag(vec)%*% mat), (vec%d*%mat), tol=1e-8 )

test.for.zero (mat2 %*% vec, mat2%d*%vec, tol=1e-8 )

test.for.zero( diag(vec)%*% vec2, vec%d*%vec2,tol=1e-8)

options(echo=TRUE)
