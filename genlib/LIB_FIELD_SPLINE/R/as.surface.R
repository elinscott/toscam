"as.surface" <-
function (obj, z, order.variables = "xy") 
{
#
    if (is.list(obj)){
      grid.list<- obj}

    if( is.matrix( obj)){
       grid.list<- attr( obj, "grid.list")}
#
#  OK now have a grid, parse this to figure 
#  nx and ny the x and y sequences and extract names
#
    hold<- parse.grid.list( grid.list, order.variables= "xy")          

#
# note that coercing z to a matrix is just reformatting 
# using the standard ordering. 
#
# output list is all the grid stuff and the matrix z. 
 
    c( hold, list(z = matrix(z, ncol = hold$ny, nrow = hold$nx)) )
}

