"unscale" <-
function (x, x.center, x.scale) 
{
    x <- scale(x, center = FALSE, scale = 1/x.scale)
    x <- scale(x, center = -x.center, scale = FALSE)
    x
}
