"radbas.constant" <-
function (m, d) 
{

# local gamma function to avoid imprecision warnings for negative arguments. 
gamma.local <- function(x){
 if( x<0){
      temp<- 1
      while( x<0){
        temp<- temp*x
        x<- x+1
      }
 return(gamma(x)/temp)}
else{
gamma(x)}
}


    if (d%%2 == 0) {
        Amd <- (((-1)^(1 + m + d/2)) * (2^(1 - 2 * m)) * 
           (pi^(-d/2)))/(gamma(m) * gamma.local(m - d/2 + 1))
    }
    else {
        Amd <- (gamma.local(d/2 - m) * (2^(-2 * m)) * (pi^(-d/2)))/gamma(m)
    }
    Amd
}
