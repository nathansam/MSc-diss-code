#' CosinorAnalysis:
#' @description Fits cosinor models to transcriptomics data and plots the
#'  best-fitting models using ggplot2.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param period The period of rhythmicity which is being tested for.
#'  Defaults to 24 (circadian).
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param threshold Set significance value for which genes are considered
#'  significant and thus plotted. Defaults to 6e-07.
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer.
#'  Defaults to TRUE.
#' @param df Logical. If TRUE a dataframe containing the results of the cosinor
#'  analysis will be returned. Defaults to TRUE.
#' @param adj String. p-value adjustment. Defaults to 'bonferonni'. 'none' is
#'  also supported.
#' @param progress Logical. If TRUE then a loading bar will be printed to denote
#'  progress.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe
#'  containing gene name and p values from F-test ranking of cosinor models.
#' @examples
#' cosinor_results <- CosinorAnalysis(Laurasmappings)
#'
#' @export


CosinorAnalysis <- function(dataset, period = 24, timelag = 6, threshold = 0.05,
                            adj = "bonferroni", progress = TRUE, save = FALSE, print = TRUE,
                            df = TRUE) {
    expected_adj <- c("bonferroni", "Bonferroni", "none")

    if (adj %in% expected_adj == FALSE) {
        stop(paste("The adjustment method ", adj, " is not recognized"))
    }


    genenumber <- nrow(dataset)  # Number of genes in the dataset
    if (adj == "bonferroni" | adj == "Bonferroni") {
        threshold <- threshold/genenumber # Calculate threshold for 'bonferroni'
    }
    pvalues <- rep(0, genenumber)  # Init list of p-values
    # First column : gene name, second column : pvalue
    cosinor.pvalue.df <- data.frame(sample = dplyr::select(dataset, 1),
                                      pVal = pvalues)
    # Calculate vector of time values
    timevector <- CircadianTools::MakeTimevector(dataset)
    loading_values <- CircadianTools::LoadingGen(genenumber)

    for (i in 1:genenumber) {
        if (progress == TRUE){
            CircadianTools::LoadingPrint(i, loading_values) # Show progress
        }

        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)

        genename <- genematrix[1, 1] # Get the genename
        genematrix <- genematrix[-1] # Get the activity readings
        genematrix <- t(genematrix) # Transpose so each row is a measurement
        # Apply time lag and create datatframe of time and activity
        geneexpression <- data.frame(timevector - timelag, genematrix)
        names(geneexpression) <- c("timevector", "activity") # Column names
        # Fit cosinor model
        cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector),
                                    period = period, data = geneexpression)
        # Get p-value of F-test
        cosinor.pvalue.df[i, 2] <- cosinor2::cosinor.detect(cosinormodel)[4]
        if (cosinor2::cosinor.detect(cosinormodel)[4] < threshold) {
            if (adj == "bonferroni" | adj == "Bonferroni") {
                plot_title <- paste("Gene=", genename, ", P-value=",
                              round(cosinor2::cosinor.detect(cosinormodel)[4] *
                  genenumber, 10))
            }

            if (adj == "none") {
                plot_title <- paste("Gene=", genename, ", P-value=",
                              round(cosinor2::cosinor.detect(cosinormodel)[4],
                  10))
            }

            p <- CircadianTools::ggplot.cosinor.lm(cosinormodel,
                                endtime = tail(timevector, n = 1) -timelag)
            p <- p + ggplot2::geom_point(ggplot2::aes(y = activity,
                    x = timevector), data = geneexpression, size = 3,
                    alpha = 0.5, color = "#39A5AE")
            p <- p + ggplot2::ggtitle(plot_title) + ggplot2::theme_bw()
            p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
            p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))
            p <- p + ggplot2::xlab("Time (hours)")
            p <- p + ggplot2::ylab("Trancripts Per Million (TPM)")


            if (save == TRUE) {
                ggplot2::ggsave(paste("Cosinor_", genename, ".png"), p,
                                width = 10, height = 4.5, units = "in")
            }
            if (print == TRUE) {
                print(p)
            }
        }
    }
    if (df == TRUE) {
        return(cosinor.pvalue.df)
    }
}


#' CosinorAnalysisPar:
#' @description Parallel Implementation of \link{CosinorAnalysis}. Fits cosinor
#'  models to transcriptomics data and plots the best-fitting models using
#'  ggplot2.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param period The period of rhythmicity which is being tested for. Defaults
#'  to 24 (circadian).
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param threshold Set significance value for which genes are considered
#'  significant and thus plotted. Defaults to 6e-07
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer.
#'  Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the results of the cosinor
#'  analysis will be returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe
#'  containing gene name and p values from F-test ranking of cosinor models
#' @examples
#' cosinor.results <- CosinorAnalysisPar(Laurasmappings)
#'
#' @export


CosinorAnalysisPar <- function(dataset, period = 24, nthreads = NULL,
                               timelag = 6) {
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }
    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }


    dataset <- CircadianTools::GeneClean(dataset)
    genenumber <- nrow(dataset)  # Number of genes in the dataset
    pvalues <- rep(0, genenumber)  # Init list of p-values
    # First column gene name, second:pvalue
    cosinor.pvalue.df <- data.frame(sample = dplyr::select(dataset, 1),
                                    pVal = pvalues)
    timevector <- CircadianTools::MakeTimevector(dataset)
    loading_values <- CircadianTools::LoadingGen(genenumber)


    cl <- parallel::makeForkCluster(nthreads)
    doParallel::registerDoParallel(cl)
    cosinor.pvalue.df <- foreach::foreach(i = 1:genenumber,
                                          .combine = rbind) %dopar% {
        CircadianTools::LoadingPrint(i, loading_values)

        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        sample <- genematrix[1, 1]
        genematrix <- genematrix[-1]
        genematrix <- t(genematrix)
        geneexpression <- data.frame(timevector - timelag, genematrix)
        names(geneexpression) <- c("timevector", "activity")
        cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector),
                                        period = period, data = geneexpression)
        pVal <- cosinor2::cosinor.detect(cosinormodel)[4]
        data.frame(sample, pVal)
    }

    return(cosinor.pvalue.df)
}


#' CosinorPlot:
#' @description Fits a cosinor model to a given gene in a given dataset and
#'  plots the model
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param period The period of rhythmicity which is being tested for. Defaults
#'  to 24 (circadian).
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#'  TRUE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' CosinorPlot('comp100000_c0_seq2', Laurasmappings )
#'
#' @export


CosinorPlot <- function(genename, dataset, timelag = 6, period = 24,
                        print = TRUE, save = FALSE) {
    genematrix <- subset(dataset, dataset[1] == genename)
    timevector <- CircadianTools::MakeTimevector(genematrix)  # Makes timevector

    genematrix <- genematrix[-1]  # Remove gene name from subset
    genematrix <- t(genematrix)
    geneexpression <- data.frame(timevector - timelag, genematrix)
    names(geneexpression) <- c("timevector", "activity")
    cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector),
                                     period = period, data = geneexpression)
    p <- CircadianTools::ggplot.cosinor.lm(cosinormodel, endtime = 21)
    p <- p + ggplot2::geom_point(ggplot2::aes(y = activity, x = timevector),
            data = geneexpression, size = 3, alpha = 0.5, color = "#008dd5")
    p <- p + ggplot2::ggtitle(paste("Gene=", genename, ", P-value=",
                     round(cosinor2::cosinor.detect(cosinormodel)[4], 10)))
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))
    p <- p + ggplot2::xlab("Time (hours)")
    p <- p + ggplot2::ylab("Trancripts Per Million (TPM)")

    if (save == TRUE) {
        ggplot2::ggsave(paste("Cosinor_", genename, ".png"), p, width = 10,
                        height = 4.5, units = "in")
    }
    if (print == TRUE) {
        return(p)
    }
}


#' CosinorSignificantPlot :
#' @description Prints or saves the genes found to be most significant by
#'  \link{CosinorAnalysis}.
#' @param results A dataframe generated by \code{CosinorAnalysis}.
#' @param dataset A transcriptomics dataset which was used with
#'  \code{CosinorAnalysis} to generate the results dataframe. First columns
#'   should be gene names. All other columns should be expression levels.
#' @param number The number of most significant genes printed or saved.
#' @param period The period of rhythmicity which is being tested for. Defaults
#'  to 24 (circadian).
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer.
#'  Defaults to TRUE.
#' @param path The path to be used for saving plots to. Only used if save is
#'  TRUE. If not specified then path will be based off name of the
#'  \code{results} object.
#' @return Prints or saves ggplot2 object(s).
#' @examples
#' cosinor.results <- CosinorAnalysis(Laurasmappings)
#' CosinorSignificantPlot(cosinor.results, Laurasmappings, number = 15, period = 24 ,save = TRUE)
#'
#' @export
CosinorSignificantPlot <- function(results, dataset, number = 10, period = 24,
                                   print = TRUE, save = FALSE, path = NULL) {
    # Order by most significant p-value
    results <- results[order(results$pVal), ]

    if (is.null(path) == TRUE) {
        # If path is not given then use name of results object
        path <- deparse(substitute(results))
    }

    if (save == TRUE) {
        if (dir.exists(path) == FALSE) {
   # If save==TRUE create directory for saved plots if doesn't already exist
            dir.create(path)
        }
    }

    for (i in 1:number) {
        p <- CircadianTools::CosinorPlot(as.character(results[i, 1]),
                                dataset, period = period)
        p <- p + ggplot2::ggtitle(paste("Gene = ",
                                as.character(results[i, 1]), " P-Value = ",
            as.character(results[i, 2])))

        if (print == TRUE) {
            print(p)
        }

        if (save == TRUE) {
            ggplot2::ggsave(paste("rank=", i, "Cosinor_",
                            as.character(results[i, 1]), ".png"), p,
                            path = path, width = 10, height = 4.5, units = "in")
        }
    }
}

#' CosinorResidualPlot
#' @description Fits a cosinor model to a gene and plots the residuals.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#' All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Shifts the plot to earlier in time.
#' @param period The period of rhythmicity which is being tested for.
#' Defaults to 24 (circadian).
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to
#' FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#' TRUE
#' @param path The directory to be used for saving plots to. Uses the working
#' directory by default. Not used if save=FALSE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' CosinorResidualPlot('comp99801_c1_seq1', Laurasmappings)
#'
#' @export
CosinorResidualPlot <- function(genename, dataset, timelag = 6, period = 24,
                                print = TRUE, save = FALSE,
                                path = NULL) {

    if (save == TRUE) {
        if (is.null(path) == FALSE) {
            if (dir.exists(path) == FALSE) {
                # If save==TRUE then create directory for saved plots if needed
                dir.create(path)
            }
        }
    }

    genematrix <- subset(dataset, dataset[1] == genename)
    timevector <- CircadianTools::MakeTimevector(genematrix)  # Makes timevector

    genematrix <- genematrix[-1]  # Remove gene name from subsetted data
    genematrix <- t(genematrix)
    geneexpression <- data.frame(timevector - timelag, genematrix)
    names(geneexpression) <- c("timevector", "activity")
    cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector),
                                        period = period, data = geneexpression)

    test <- data.frame(geneexpression, resids = cosinormodel$fit$residuals)

    ymax <- max(abs(test$resids)) # Find maximum
    ymax <- ceiling(ymax)

    p <- ggplot2::ggplot(ggplot2::aes(x = timevector, y = resids), data = test)
    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
                                 color = "#ffa630", size = 1.1, alpha = 1)
    p <- p + ggplot2::geom_point(size = 3.5, alpha = 0.5, color = "#008dd5")
    p <- p + ggplot2::xlab("Time (Hours)")
    p <- p + ggplot2::ylab("Residuals") + ggplot2::theme_bw()
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))
    p <- p + ggplot2::ggtitle(paste("Residual Plot for ", genename, " with ",
                                    period, " Hour Period", sep = ""))
    p <- p + ggplot2::ylim(-1 * ymax, ymax)

    if (save == TRUE) {
        ggplot2::ggsave(paste(genename, ".png", sep = ""), p, path = path,
                        width = 10, height = 4.5, units = "in")
    }

    if (print == TRUE) {
        return(p)
    }
}


#' CosinorResidualDatasetPlot
#' @description Fits a cosinor model and plots the residuals for multiple genes
#' in a dataset.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#' All other columns should be expression levels.
#' @param timelag Shifts the plot to earlier in time.
#' @param period The period of rhythmicity which is being tested for.
#' Defaults to 24 (circadian).
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to
#' FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#' TRUE
#' @param path The directory to be used for saving plots to. Uses the working
#' directory by default. Not used if save=FALSE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' CosinorResidualDatasetPlot(Laurasmappings[1 : 5, ])
#'
#' @export
CosinorResidualDatasetPlot <- function(dataset, timelag = 6, period = 24,
                                       print = TRUE, save = FALSE, path = NULL){

    genes <- dataset[,1] # Find gene names

    for (i in genes){
        CosinorResidualPlot(i, dataset = dataset, timelag = timelag,
                       period = period, print = print, save = save, path = path)
    }
}

