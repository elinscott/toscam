"make.surface.grid" <-
function (grid.list)
{
#
# the old fields version of make.surface.grid was complicated
# and we believe rarely used. 
# this current function 
# is essentially a single line replacement 
# 
# but adds an attribute for the grid matrix to carry 
# and carries along the names of the grid.list variables.
# along the information as to how it was created. 
# see as.surface 

  temp<-as.matrix( expand.grid( grid.list))

# wipe out row names
  dimnames(temp)<-list( NULL, names(grid.list))

# set attribute                 
  attr(temp, "grid.list")<- grid.list

temp
}

