#' @export
LoadingGen <- function(genenumber) {
    loading_values <- rep(0, 20)
    for (i in 1:20) {
        loading_values[i] <- round(genenumber/20 * i)
    }
    return(loading_values)
}

#' @export
LoadingPrint <- function(iteration, loading_values) {
    if (iteration == 1) {
        cat(crayon::red(noquote("This may take a while if using a large dataset! \n")))
        cat(noquote("Progress: \n"))
        cat(noquote(paste(replicate(20, "□"), collapse = "")))
        cat(noquote("\n"))
    }
    if (iteration %in% loading_values) {
        position <- match(iteration, loading_values)
        hashes <- paste(replicate(position, "■"), collapse = "")
        dashes <- paste(replicate(20 - position, "□"), collapse = "")
        cat(noquote(paste("\r", hashes, dashes, "\n", sep = "")))
        flush.console()
    }
}
