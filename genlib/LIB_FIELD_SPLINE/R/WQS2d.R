"WQS2d" <-
function (x, transpose = FALSE) 
{
    if (!transpose) 
        t(WQS(t(WQS(x))))
    else {
        WQS.T(t(WQS.T(t(x))))
    }
}
