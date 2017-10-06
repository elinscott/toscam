"WQS.basis" <-
function (N, cut.n = 8) 
{
Wtransform( diag(1,N), inv=TRUE,cut.min=cut.n)
}
