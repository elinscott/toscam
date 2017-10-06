"Exp.earth.cov" <-
function (x1, x2, theta = 1) 
{
    exp(-rdist.earth(x1, x2)/theta)
}
