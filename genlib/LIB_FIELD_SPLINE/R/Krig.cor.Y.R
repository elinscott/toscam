Krig.cor.Y<- function(obj, verbose=FALSE ){

# subtract mean
if( !is.na( obj$mean.obj[1])){
Y<- obj$y - predict(obj$mean.obj, obj$x)
}

# divide by sd
if( !is.na( obj$sd.obj[1])){
Y<- Y/predict(obj$sd.obj,obj$x)
}

Y
 

}

