#' CorAnalysis:
#' @description Ranks correlation between a given gene and all other genes in a
#'  dataset. Plots both the given gene and highly correlated genes for a given
#'  correlation value
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param genename the name of a gene intended for comparison with all other
#' genes in the dataset. Must be a string.
#' @param threshold Set correlation threshold value for which genes are
#'  considered significant and thus plotted. Defaults to 0.9
#' @param lag Setting any value other than 0 allows a gene to be correlated with
#'  lagged genes in the dataset. The number denotes the number of timesteps to
#'  lag by.
#' @param average The average to be used for comparing the time points. Either
#'  'median' or 'mean'.
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders highly correlated genes in the plot
#' viewer. Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the correlations of the
#'  given gene with all genes in the dataset is returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe
#'  containing gene names and correlation values
#' @examples
#' cor_results <- coranalysis('comp100002_c0_seq2', Laurasmappings)
#'
#' @export


CorAnalysis <- function(genename, dataset, threshold = 0.9, average = "median",
                        lag = 0, save = FALSE, print = TRUE, df = TRUE) {

    if (save == TRUE) {
        directory <- paste("cor_", genename)
        if (dir.exists(directory) == FALSE) {
            dir.create(directory)
        }
    }

    genenumber <- nrow(dataset)  # Number of genes in the dataset
    # First column gene name, second column correlation value
    cor.df <- data.frame(sample = dplyr::select(dataset, 1),
                         corvalues = rep(0, genenumber))
    # Create vector of time values
    timevector <- CircadianTools::MakeTimevector(dataset)
    # Used for the loading bar
    loading_values <- CircadianTools::LoadingGen(genenumber)
    selectedgene <- CircadianTools::ActivitySelect(genename, dataset)
    selectedgenedf <- data.frame(timevector, selectedgene)
    names(selectedgenedf) <- c("timevector", "activity")


    selectedaverage.list <- rep(0, length((unique(timevector))))
    count <- 1
    for (i in unique(timevector)) {
        genesub <- subset(selectedgenedf, timevector == i, select = activity)
        if (average == "mean") {
            selectedaverage.list[[count]] <- (mean(genesub$activity))
        }
        if (average == "median") {
            selectedaverage.list[[count]] <- (median(genesub$activity))
        }
        count = count + 1
    }

    if (lag > 0) {
        selectedaverage.list <- tail(selectedaverage.list,
                                     n = length(selectedaverage.list) - lag)
    }

    if (lag < 0) {
        selectedaverage.list <- head(selectedaverage.list,
                                     n = length(selectedaverage.list) - lag)
    }

    for (i in 1:genenumber) {
        CircadianTools::LoadingPrint(i, loading_values)

        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        compgenename <- paste(dataset[i, 1])
        genematrix <- CircadianTools::ActivitySelect(i, dataset)



        selectedgenedf <- data.frame(timevector, genematrix)

        names(selectedgenedf) <- c("timevector", "activity")


        compaverage.list <- rep(0, length((unique(timevector))))
        count <- 1
        for (j in unique(timevector)) {
            compgenesub <- subset(selectedgenedf, timevector == j,
                                  select = activity)
            if (average == "mean") {
                compaverage.list[[count]] <- (mean(compgenesub$activity))
            }
            if (average == "median") {
                compaverage.list[[count]] <- (median(compgenesub$activity))
            }
            count = count + 1
        }

        if (lag > 0) {
            compaverage.list <- head(compaverage.list,
                                     n = length(compaverage.list) - lag)
        }
        if (lag < 0) {
            compaverage.list <- tail(compaverage.list,
                                     n = length(compaverage.list) - lag)
        }


        correlation <- cor(selectedaverage.list, compaverage.list)
        cor.df[i, 2] <- correlation
        if (save == TRUE || print == TRUE) {
            if (correlation > threshold) {
                if (correlation != 1) {


                  myplot <- compplot(as.character(genename), compgenename,
                                     dataset)
                  myplot <- myplot + ggplot2::ggtitle(paste("Correlation = ",
                                                            correlation))

                  if (save == TRUE) {
                    ggplot2::ggsave(paste("Cor_", genename, "_",
                                          compgenename, ".png"), myplot,
                                    path = directory, width = 10, height = 4.5,
                                    units = "in")
                  }
                  if (print == TRUE) {
                    print(myplot)
                  }
                }
            }
        }
    }
    if (df == TRUE) {
        return(cor.df)
    }
}

#' CorAnalysisCluster
#' @description Correlates the average activity of a cluster with the average
#'  activity of every other cluster.
#' @param cluster.no The ID for the cluster which will be compared with all
#'  other clusters in the dataset
#' @param cluster.dataset A transcriptomics dataset where the final column
#'  details the cluster the gene belongs to. First column should be gene names.
#'  All remaining columns should be expression levels.
#' @param lag Setting any value other than 0 allows a cluster to be correlated
#'  with lagged clusters in the dataset. The number denotes the number of
#'  timesteps to lag by.
#' @param nthreads The number of threads to be used for parallel computations.
#'  Defaults to the maximum number of threads available.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, 50)
#' cor.df <- CorAnalysisCluster(3, pam.df)
#'
#' @export

CorAnalysisCluster <- function(cluster.no, cluster.dataset, lag = 0,
                               nthreads = NULL) {


    # Get time profile for the cluster specified
    main.time.profile <- CircadianTools::ClusterTimeProfile(cluster.no,
                                cluster.dataset, nthreads = nthreads)

    if (lag > 0) {
        main.time.profile <- tail(main.time.profile,
                     n = length(main.time.profile) - lag)  # Lag if required
    }
    if (lag < 0) {
        main.time.profile <- head(main.time.profile,
                    n = length(main.time.profile) - lag)  # Lag if required
    }

    cluster.quantity <- max(cluster.dataset$cluster)  # Number of clusters
    correlation.df <- data.frame(seq(1, cluster.quantity),
                                     rep(0, cluster.quantity))
    colnames(correlation.df) <- c("cluster", "correlation")

    for (j in 1:cluster.quantity) {
        # Get time profile of cluster being compared with the main cluster
        comp.time.profile <- CircadianTools::ClusterTimeProfile(j,
                            cluster.dataset, nthreads = nthreads)

        if (lag > 0) {
            comp.time.profile <- head(comp.time.profile,
                        n = length(comp.time.profile) - lag)  # Lag if required
        }
        if (lag < 0) {
            comp.time.profile <- tail(comp.time.profile,
                        n = length(comp.time.profile) - lag)  # Lag if required
        }

        # Calculate correlation
        compcor <- cor(main.time.profile, comp.time.profile)
        # Add correlation to dataframe
        correlation.df[j, 2] <- compcor
    }
    return(correlation.df)
}

#' CorAnalysisDataset:
#'
#' @description Correlates every gene in a dataset with every other gene in the
#' same dataset. Allows a timelag between genes to be correlated.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param average The average to be used when comparing the time points. Either
#'  'median' or 'mean'.
#' @param lag Setting any value other than 0 allows a gene to be correlated with
#'  lagged genes in the dataset. The number denotes the number of timesteps to
#'  lag by.
#' @param nthreads The number of threads to be used for parallel computations.
#'  Defaults to the maximum number of threads available.
#' @param save Logical. If TRUE, the dataframe of correlations for the dataset
#'  is saved as a .csv file.
#' @param filename filename for saved csv file. Only used if save=TRUE. If not
#'  specified then the dataset object name is used.
#' @return A dataframe of correlation values. The column genes represent the
#'  original genes whilst the rows represent lagged genes.
#' @examples
#' subdf<-TFilter(Laurasmappings)
#' cordf <- CorAnalysisDataset(subdf, lag=1,filename='cor_tfiltered')
#'
#' @export
CorAnalysisDataset <- function(dataset, average = "median", lag = 0,
                               nthreads = NULL, save = TRUE, filename = NULL) {
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    if (is.null(filename) == TRUE) {
   # If a filename isn't specified then the name of the dataframe object is used
        filename <- deparse(substitute(dataset))
    }

    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    genenames <- as.vector(dataset[, 1])  # Vector of names for every gene

    correlationdf <- foreach::foreach(i = 1:length(genenames),
                                      .combine = cbind) %dopar% {

     # Calculate correlation for the ith gene with all ( possibly lagged) genes.
        results <- CircadianTools::CorAnalysis(genename = genenames[i],
                    dataset = dataset, lag = lag, average = average,
                      print = FALSE)

        # Add the list of correlations as a column to the correlation dataframe
        data.frame(results[, 2])
    }

    rownames(correlationdf) <- genenames  # Give the columns the genenames
    colnames(correlationdf) <- genenames  # Give the rows the genenames
    if (save == TRUE) {
        # Write as a csv file if save=TRUE
        write.csv(correlationdf, paste(filename, ".csv", sep = ""))
    }
    parallel::stopCluster(cl)
    return(correlationdf)  # Return the correlation dataframe.
}

#' CorAnalysisClusterDataset:
#'
#' @description Correlates the average activity of each cluster with every other
#'  cluster in a dataset.
#' @param cluster.dataset A transcriptomics dataset where the final column
#' details the cluster the gene belongs to. First column should be gene names.
#' All remaining columns should be expression levels.
#' @param lag Setting any value other than 0 allows a gene to be correlated with
#'  lagged genes in the dataset. The number denotes the number of timesteps to
#'  lag by.
#' @param nthreads The number of threads to be used for parallel computations.
#'  Defaults to the maximum number of threads available.
#' @param save Logical. If TRUE, the dataframe of correlations for the dataset
#'  is saved as a .csv file.
#' @param filename filename for saved csv file. Only used if save=TRUE. If not
#'  specified then the dataset object name is used.
#' @return A dataframe of correlation values. The column genes represent the
#'  original clusters whilst the rows represent lagged clusters.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filterdf)
#' cor.df <- CorAnalysisClusterDataset(pam.df)
#'
#' @export

CorAnalysisClusterDataset <- function(cluster.dataset, lag = 0, nthreads = NULL,
                                      save = TRUE, filename = NULL) {

    if (is.null(filename) == TRUE) {
   # If a filename isn't specified then the name of the dataframe object is used
        filename <- deparse(substitute(cluster.dataset))
    }

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)


    clusters <- unique(cluster.dataset$cluster)  # Vector of cluster numbers

    correlationdf <- foreach::foreach(i = 1:length(clusters),
                                      .combine = cbind) %dopar% {
        temp.df <- CircadianTools::CorAnalysisCluster(i, cluster.dataset,
                                                      lag = lag, nthreads = 1)
        temp.df[, 2]
    }

    rownames(correlationdf) <- clusters  # Give the columns the genenames
    colnames(correlationdf) <- clusters  # Give the rows the genenames

    if (save == TRUE) {
        # Write as a csv file if save=TRUE
        write.csv(correlationdf, paste(filename, ".csv", sep = ""))
    }
    parallel::stopCluster(cl)

    return(correlationdf)
}

#' CorAnalysisPar:
#' @description Parallel Implementation of \link{CorAnalysis}. Ranks correlation
#'  between a given gene and all over genes in a dataset. Plots both the given
#'  gene and highly correlated genes for a given correlation value
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param genename the name of a gene intended for comparison with all other
#'  genes in the dataset. Must be a string.
#' @param lag Setting any value other than 0 allows a gene to be correlated with
#'  lagged genes in the dataset. The number denotes the number of timesteps
#'  to lag by.
#' @param average The average to be used for comparing the time points. Either
#'  'median' or 'mean'.
#' @param nthreads Number of processor threads for the process. If not specifed
#'  then the maximum number of logical cores are used.
#' @return Returns dataframe containing gene names and correlation values
#' @examples
#' cor_results <- CorAnalysisPar('comp100002_c0_seq2', Laurasmappings)
#'
#' @export


CorAnalysisPar <- function(genename, dataset, lag = 0, average = "median",
                           nthreads = NULL) {

    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }


    genenumber <- nrow(dataset)  # Number of genes in the dataset
    # First column gene name, second column correlation value
    cor.df <- data.frame(sample = dplyr::select(dataset, 1),
                         corvalues = rep(0, genenumber))
    # Create vector of time values
    timevector <- CircadianTools::MakeTimevector(dataset)
    selectedgene <- as.vector(CircadianTools::ActivitySelect(genename, dataset))
    selectedgenedf <- data.frame(timevector, selectedgene)
    colnames(selectedgenedf) <- c("timevector", "activity")


    selectedaverage.list <- rep(0, length((unique(timevector))))
    count <- 1
    for (i in unique(timevector)) {
        genesub <- subset(selectedgenedf, timevector == i, select = activity)
        if (average == "mean") {
            selectedaverage.list[[count]] <- (mean(genesub$activity))
        }
        if (average == "median") {
            selectedaverage.list[[count]] <- (median(genesub$activity))
        }
        count = count + 1
    }

    if (lag > 0) {
        selectedaverage.list <- tail(selectedaverage.list,
                                     n = length(selectedaverage.list) - lag)
    }

    if (lag < 0) {
        selectedaverage.list <- head(selectedaverage.list,
                                     n = length(selectedaverage.list) - lag)
    }


    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    cor.df <- foreach::foreach(i = 1:genenumber, .combine = rbind) %dopar% {

        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        compgenename <- paste(dataset[i, 1])
        genematrix <- CircadianTools::ActivitySelect(i, dataset)

        selectedgenedf <- data.frame(timevector, genematrix)

        names(selectedgenedf) <- c("timevector", "activity")


        compaverage.list <- rep(0, length((unique(timevector))))
        count <- 1
        for (j in unique(timevector)) {
            compgenesub <- subset(selectedgenedf, timevector == j,
                                  select = activity)
            if (average == "mean") {
                compaverage.list[[count]] <- (mean(compgenesub$activity))
            }
            if (average == "median") {
                compaverage.list[[count]] <- (median(compgenesub$activity))
            }
            count = count + 1
        }
        if (lag > 0) {
            compaverage.list <- head(compaverage.list,
                                     n = length(compaverage.list) - lag)
        }
        if (lag < 0) {
            compaverage.list <- tail(compaverage.list,
                                     n = length(compaverage.list) - lag)
        }

        correlation <- cor(selectedaverage.list, compaverage.list)
        data.frame(compgenename, correlation)

    }

    parallel::stopCluster(cl)
    return(cor.df)
}

#' CorSignificantPlot :
#' @description Prints or saves the genes found to be most significant by
#'  \link{CorSignificantPlot} or \link{CorAnalysisPar}.
#'
#' @param results A dataframe generated by \code{CorAnalysis}.
#' @param dataset A transcriptomics dataset which was used with
#' \code{CorAnalysis} to generate the results dataframe. First columns should be
#'  gene names. All other columns should be expression levels.
#' @param number The number of most significant genes printed or saved.
#' @param period The period of rhythmicity which is being tested for. Defaults
#'  to 24 (Circadian).
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to
#'  FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer.
#'  Defaults to TRUE
#' @return Prints or saves ggplot2 object(s).
#' @examples
#' cor_results <- CorAnalysis('comp100002_c0_seq2',Laurasmappings)
#' CorSignificantPlot(cor_results, Laurasmappings, number = 15, save=TRUE, negative=FALSE)
#'
#' @export
CorSignificantPlot <- function(results, dataset, number = 10, print = TRUE,
                               save = FALSE, negative = FALSE) {
    # Order by most positive correlation
    results <- results[order(results$correlation, decreasing = TRUE), ]
    gene1 <- as.character(results[1, 1])
    if (save == TRUE) {
        directory <- paste("cor_", gene1)
        if (dir.exists(directory) == FALSE) {
            dir.create(directory)
        }
    }
    if (negative == FALSE) {
        for (i in 2:(number + 1)) {
            myplot <- CircadianTools::CompPlot(as.character(gene1),
                                        as.character(results[i, 1]), dataset)
            myplot <- myplot + ggplot2::ggtitle(paste(" Correlation = ",
                                                  as.character(results[i, 2])))

            if (print == TRUE) {
                print(myplot)
            }

            if (save == TRUE) {
                ggplot2::ggsave(paste("rank=", i - 1, "Cor_",
                    as.character(results[i, 1]), ".png"), myplot,
                       path = directory, width = 10, height = 4.5, units = "in")
            }
        }
    }
    if (negative == TRUE) {
        results <- results[order(results$correlation, decreasing = FALSE), ]
        for (i in 1:number) {
            myplot <- CircadianTools::CompPlot(as.character(gene1),
                         as.character(results[i, 1]), dataset)
            myplot <- myplot + ggplot2::ggtitle(paste("Gene = ",
                            as.character(results[i, 1]), " Cor = ",
                                   as.character(results[i, 2])))

            if (print == TRUE) {
                print(myplot)
            }

            if (save == TRUE) {
                ggplot2::ggsave(paste("rank=", i, "Cor_",
                     as.character(results[i, 1]), ".png"), myplot,
                       path = directory, width = 10, height = 4.5, units = "in")
            }
        }
    }
}

