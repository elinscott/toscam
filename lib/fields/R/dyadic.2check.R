dyadic.2check<- function( m,n,cut.p=2){
# checks that n is of the form 
# n=p*2^m where p <= cut.p
m2<- as.integer(m)
n2<- as.integer(n)

while( (n2 > cut.p) &(m2 > cut.p)) {
if((m2%%2!=0)|(n2%%2!=0) ) {
cat(n,"and" , m, 
"must equal p*2^L where p is less than or equal to ", cut.p
, fill=TRUE)
return(FALSE)
}
m2<- m2/2
n2<- n2/2

}
return(TRUE)
}

  
