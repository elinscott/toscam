"fast.1way" <-
function (lev, y, w = rep(1, length(y))) 
{
    N <- length(y)
# ordered unique values of lev
    tags <- lev[!duplicated(lev)]
# lev are now integer tags
    lev <- match(lev, tags)

# add together weights with same lev
    w.means<- c(tapply( w, lev, sum))

# find weighted means for each lev
    means <- c(tapply( y*w, lev, sum)/w.means)

# find SS
    SSE <- sum(w * (y - means[lev])^2)
    MSE <- SSE/(length(y) - length(means))

  list(n=length(means) , means = means, SSE = SSE, w.means = w.means, 
        MSE = MSE, lev = lev, tags=tags)
}

