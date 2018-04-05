"WQS2di" <-
function (x, transpose = FALSE) 
{
    if (!transpose) 
        t(WQSi(t(WQSi(x))))
    else {
        WQSi.T(t(WQSi.T(t(x))))
    }
}
