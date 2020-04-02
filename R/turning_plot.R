#' TurningPlot:
#' @description  Fits a spline to a given gene in a given dataset. Finds the
#'  turning points. Plots the turning points and spline.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Lags the time. Usually desired if wanting to start from t = 0.
#'  Defaults to 0.
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#'  TRUE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' TurningPlot('comp100000_c0_seq2',LauraSingleMap)
#'
#' @export


TurningPlot <- function(genename, dataset, timelag = 0, print = TRUE,
                        save = FALSE) {
    genematrix <- subset(dataset, dataset[1] == genename)
    # Make vector of time values
    timevector <- CircadianTools::MakeTimevector(genematrix) - timelag
    genematrix <- genematrix[-1]  # Remove gene name from subset

    genesplinefunc <- splinefun(timevector, genematrix)
    x <- seq(timevector[1], tail(timevector, n = 1), by = 0.001)
    y <- genesplinefunc(x, deriv = 0)  # Deriv=0 gives the spline itself
    z <- genesplinefunc(x, deriv = 1)  # First derivative of the spline

    spline.df <- data.frame(x, y, z)

    # This outputs a vector of the turning points
    # i.e. the points where the derivative=0
    turning.points <- rootSolve::uniroot.all(genesplinefunc,
                                             interval = c(6, 27), deriv = 1)

    p <- ggplot2::ggplot(ggplot2::aes(x = x, y = y), data = spline.df)
    p <- p + ggplot2::geom_line(color = "#008dd5", size = 2)
    p <- p + ggplot2::theme_bw() + ggplot2::xlab("Time(hours)")
    p <- p + ggplot2::ylab("Transcripts Per Million (TPM)")
    p <- p + ggplot2::ggtitle(paste("Gene=", genename, ", Mean difference=",
                                    round(mean(diff(turning.points)),2),
                                    " SD=",
                                    round(sd(diff(turning.points)), 2)))
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))

    for (i in 1:length(turning.points)) {
        p <- p + ggplot2::geom_vline(xintercept = turning.points[i], size = 1.2,
                                     alpha = 0.75, color = "#ba1200")
    }

    p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))


    if (save == TRUE) {
        ggplot2::ggsave(paste("Turning_", genename, ".png"), p, width = 10,
                        height = 4.5, units = "in")
    }

    if (print == TRUE) {
        return(p)
    }
}
