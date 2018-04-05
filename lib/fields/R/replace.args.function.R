"replace.args.function" <-
function (fun, ...) 
{
    temp <- list(...)
    ntemp <- names(temp)
    fnames <- names(fun)
    if (length(temp) > 0) {
        for (k in 1:length(ntemp)) {
            if (!is.na(match(ntemp[k], fnames))) {
                fun[ntemp[k]] <- temp[ntemp[k]]
            }
        }
    }
    as.function(fun)
}
