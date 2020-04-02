#' BasicPlot:
#' @description Plots activity data as points and average activity as lines
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Shifts the plot to earlier in time.
#' @param method How should the average activity for each time point be
#'  calculated? 'median' or 'mean'
#' @param points Logical. If FALSE then each observation will not be plotted
#' as points.
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#'  TRUE
#' @param path The directory to be used for saving plots to. Uses the working
#'  directory by default. Not used if save = FALSE
#' @return A ggplot2 object
#' @examples
#' basicplot('comp101252_c0_seq2',Laurasmappings)
#'
#' @export

BasicPlot <- function(genename, dataset, timelag = 0, method = "median",
                      points = TRUE, print = TRUE, save = FALSE, path = NULL) {
    if (save == TRUE) {
        if (is.null(path) == FALSE) {
            if (dir.exists(path) == FALSE) {
                # If save==TRUE then create directory for saved plots if needed
                dir.create(path)
            }
        }
    }

    genematrix <- subset(dataset, dataset[1] == genename)  # Select gene
    # Makes vector of time values
    timevector <- (CircadianTools::MakeTimevector(genematrix) - timelag)
    genematrix <- t(genematrix[-1])  # Remove gene name

    genedata <- data.frame(timevector, genematrix)
    names(genedata) <- c("timevector", "activity")

  # Make a list to hold the averaged values for gene activity at each time point
    average.list <- rep(0, length((unique(timevector))))

    count <- 1
    for (i in timevector) {
        genesub <- subset(genedata, timevector == i, select = activity)
        if (method == "mean")
            {
                average.list[[count]] <- (mean(genesub$activity))
            }  #compute the mean of the time points
        if (method == "median")
            {
                average.list[[count]] <- (median(genesub$activity))
            }  #compute the median of the time points
        count = count + 1
    }

    genedata <- cbind(genedata$timevector, genedata$activity, average.list)
    genedata <- data.frame(genedata)
    names(genedata) <- c("timevector", "activity", "new.average")
    p <- ggplot2::ggplot(data = genedata,
                         ggplot2::aes(x = timevector,y = activity))
    p <- p + ggplot2::geom_line(ggplot2::aes(x = timevector, y = new.average),
                                size = 1, color = "#412d6b")

    if (points == TRUE) {
        p <- p + ggplot2::geom_point(size = 3, alpha = 0.5, color = "#008dd5")
    }
    p <- p + ggplot2::xlab("Time (Hours)")
    p <- p + ggplot2::ylab("Transcripts Per Million (TPM)")
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))
    p <- p + ggplot2::ggtitle(paste("Gene = ", genename))

    if (save == TRUE) {
        ggplot2::ggsave(paste(genename, ".png"), p, path = path, width = 10,
                        height = 4.5, units = "in")

    }

    if (print == TRUE) {
        return(p)
    }
}


#' CompPlot:
#' @description Plots two genes from a gene activity dataset
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param gene1 The name of the first gene to be plotted. Must be a string.
#' @param gene2 The name of the second gene to be plotted. Must be a string.
#' @param save Logical. If TRUE, saves plots to working directory.
#'  Defaults to FALSE.
#' @param points Logical. If FALSE then each observation will not be plotted
#'  as points.
#'
#'
#' @return Returns or saves a ggplot2 object.
#' @examples
#' CompPlot('comp100000_c0_seq3', 'comp100002_c0_seq2', Laurasmappings)
#'
#' @export

CompPlot <- function(gene1, gene2, dataset, save = FALSE, points = TRUE) {

    # Vector of time values
    timevector <- CircadianTools::MakeTimevector(dataset)

    # Get gene 1's activity values
    activity.1 <- subset(dataset, sample == gene1)
    activity.1 <- t(activity.1[-1])

    # Get gene 2's activity values
    activity.2 <- subset(dataset, sample == gene2)
    activity.2 <- t(activity.2[-1])


    combined <- cbind(timevector, activity.1, activity.2)
    combined <- data.frame(combined)
    names(combined) <- c("timevector", "gene1_activity", "gene2_activity")


    mean.list <- rep(0, length(unique(timevector)))
    count <- 1
    for (i in unique(timevector)) {
        bob <- subset(combined, timevector == i, select = gene1_activity)
        mean.list[[count]] <- (mean(bob$gene1_activity))
        count = count + 1
    }

    mean.list2 <- rep(0, length(unique(timevector)))
    count <- 1
    for (i in unique(timevector)) {
        bob <- subset(combined, timevector == i, select = gene2_activity)
        mean.list2[[count]] <- (mean(bob$gene2_activity))
        count <- count + 1
    }

    new_mean_list <- rep(mean.list, as.numeric(table(timevector)))
    new_mean_list <- data.matrix(new_mean_list)

    new_mean_list2 <- rep(mean.list2, as.numeric(table(timevector)))
    new_mean_list2 <- data.matrix(new_mean_list2)

    combined <- cbind(combined$timevector, combined$gene1_activity,
                      combined$gene2_activity, new_mean_list, new_mean_list2)

    combined <- data.frame(combined)

    names(combined) <- c("timevector", "gene1_activity", "gene2_activity",
                         "mean1", "mean2")


    p <- ggplot2::ggplot(data = combined, ggplot2::aes(timevector))
    p <- p + ggplot2::scale_x_continuous(breaks = unique(timevector))
    p <- p + ggplot2::geom_line(ggplot2::aes(y = mean1, colour = "first gene"),
                                size = 1)
    p <- p + ggplot2::xlab("Time (hours)")
    p <- p + ggplot2::ylab("Trancripts Per Million (TPM)")
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::geom_line(ggplot2::aes(y = mean2, color = "second gene"),
                                size = 1)
    if (points == TRUE) {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = gene2_activity,
                                                  colour = "second gene"),
            size = 3, alpha = 0.5)
        p <- p + ggplot2::geom_point(ggplot2::aes(y = gene1_activity,
                                                  colour = "first gene"),
            size = 3, alpha = 0.5)
    }

    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
    p <- p + ggplot2::scale_color_manual(name = "Gene",
                                        breaks = c("first gene", "second gene"),
                                        labels = c(gene1, gene2),
                                        values = c("#008dd5", "#ffa630"))

    if (save == TRUE) {
        ggplot2::ggsave(paste(gene1, gene2, "comp.png"), p, width = 10,
                        height = 4.5, units = "in")
    }

    return(p)
}


#' DatasetPlot:
#' @description Saves plots of all genes in a dataset. WARNING! Don't run on a
#'  large dataset! Intended for a filtered dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param timelag Shifts the plot to earlier in time.
#' @param points Logical. If FALSE then each observation will not be plotted
#' as points.
#' @param method How should the average activity for each time point be
#'  calculated? 'median' or 'mean'
#' @param nthreads Number of processor threads to be used for calculating the
#' distance matrix. If not specifed then the maximum number of logical cores
#' are used.
#' @param path The directory to be used for saving plots to. Uses the name of
#' the dataset object if this argument is not specified)  Not used if save=FALSE
#' @examples
#' filter.df<-CombiFilter(Laurasmappings)
#' DatasetPlot(filter.df)
#' @export

DatasetPlot <- function(dataset, timelag = 0, method = "median", points = TRUE,
                        nthreads = NULL, path = NULL) {

    if (is.null(path) == TRUE) {
   # If a filename isn't specified then the name of the dataframe object is used
        path <- deparse(substitute(dataset))
    }
    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    genenames <- dataset[, 1]
    foreach::foreach(i = 1:length(genenames)) %dopar% {
        CircadianTools::BasicPlot(genenames[i], dataset = dataset,
                                  timelag = timelag, method = method,
                                  points = points, print = FALSE, save = TRUE,
                                  path = path)
    }
}
