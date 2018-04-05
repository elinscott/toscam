"Exponential" <-
function (d , range = 1, alpha=1/range, phi=1) 
{
#
# Matern covariance function transcribed from Stein's book page 31
# nu==smoothness==.5, alpha ==  1/range
#
# GeoR parameters map to kappa==smoothness and phi == range 
# check for negative distances

if( any( d <0)) stop("distance argument must be nonnegative")

#
return( phi*exp(-d*alpha) ) 
}
