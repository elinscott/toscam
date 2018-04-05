"ceiling2" <-
function (m) 
{
    if (m < 1) 
        return(NA)
    M <- 1
    while (M < m) {
        M <- M * 2
    }
    M
}
