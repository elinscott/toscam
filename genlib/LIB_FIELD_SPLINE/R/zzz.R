".First.lib" <-
function (lib, pkg) 
{library.dynam("fields",pkg, lib)
cat("fields is loaded use help(fields) for an overview of this library ", 
"  Some name changes have been made to several common functions:",
"  exp.-> Exp. and rad. -> Rad. See help( Exp.cov)", fill=TRUE)
}
