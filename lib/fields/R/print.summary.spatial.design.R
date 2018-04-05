"print.summary.spatial.design" <-
function (x, digits = 4,...) 
{
    cat("Call:\n")
    dput(x$call)
    c1 <- "Number of design points:"
    c2 <- length(x$best.id)
    c1 <- c(c1, "Number of fixed points:")
    if (is.null(x$fixed)) 
        c2 <- c(c2, 0)
    else c2 <- c(c2, length(x$fixed))
    c1 <- c(c1, "Optimality Criterion:")
    c2 <- c(c2, round(x$opt.crit, digits))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    print(sum, quote = FALSE, digits = digits)
    other.crit <- x$other.crit
    if (length(other.crit) > 1) {
        cat("\nOptimality criteria for other designs:\n\t")
        cat(round(other.crit, digits), "\n")
    }
    cat("\nHistory:\n")
    dimnames(x$history)[[1]] <- rep("", nrow(x$history))
    print(round(x$history, digits), quote = FALSE)
    invisible(x)
}
