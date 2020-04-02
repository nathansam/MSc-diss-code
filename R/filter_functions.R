#' AnovaFilter
#' @description Filters a gene activity dataframe via ANOVA.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not
#' specified then the maximum number of logical cores are used.
#' @param threshold Set the p-value threshold for the filtering.
#' @examples
#' Laurasmappings_filtered <- AnovaFilter(Laurasmappings, nthreads=4)
#'
#' @export

AnovaFilter <- function(dataset, threshold = 0.05, nthreads = NULL) {
    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    genenumber <- nrow(dataset)
    # List of time values (repeated for replicates)
    timevector <- CircadianTools::MakeTimevector(dataset)

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)  # Register cluster

    filterdf <- foreach::foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        # Parallel for loop to create dataframe of significant genes
        # Get gene by row
        gene <- dplyr::filter(dataset, dplyr::row_number() == i)
        genematrix <- t(gene[-1])  # Remove gene name
        # Fit ANOVA model and create aov object
        tempaov <- aov(lm(as.numeric(genematrix) ~ as.factor(timevector)))
        # Get the p-value from the aov object
        pvalue <- summary(tempaov)[[1]][1, 5]
        if (pvalue < threshold) {
            # Return the gene (this gene will be combined with other significant
            # genes found in the for loop via rbind to form the dataframe)
            gene
        }
    }
    parallel::stopCluster(cl)

    rownames(filterdf) <- seq(1, nrow(filterdf))  # Rebuild the row names

    return(filterdf)
}


#' SizeFilter
#' @description Filters the genes with the smallest range from a transcriptomics
#'  dataset.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param cutoff Proportion of total number of genes (which show at least some
#'  activity) to be removed.
#' @param nthreads Number of processor threads for the filtering. If not
#'  specifeid then the maximum number of logical cores are used.
#' @examples
#' rangedf <- SizeFilter(Laurasmappings, nthreads=4)
#'
#' @export

SizeFilter <- function(dataset, cutoff = 0.1, nthreads = NULL) {
    # Remove genes with no activity
    dataset <- CircadianTools::GeneClean(dataset)
    genenumber <- nrow(dataset)
    timevector <- CircadianTools::MakeTimevector(dataset)

    # Calculate range of activity for each gene
    rangedf <- CircadianTools::GeneRange(dataset = dataset, nthreads = nthreads)
    # Order by biggest range
    rangedf <- rangedf[order(rangedf$generange, decreasing = TRUE), ]
    # Number of genes which will be returned after filtering
    rank <- round(genenumber * (1 - cutoff))
    # Selects the genes with the highest range
    rangedf <- rangedf[1:rank, ]
    # Returns genes found in rangedf with their respective activity readings
    filterdf <- CircadianTools::GeneSub(rangedf, dataset, nthreads)
    row.names(filterdf) <- seq(1, nrow(filterdf))

    return(filterdf)
}

#' ZeroFilter
#' @description Filters a transcriptomics dataset such that there is a minimum
#'  number of non-zero activity readings for each gene in the resultant dataset.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param non_zero_num The minimum number of non-zero activity readings for
#'  genes in the filtered dataset.
#' @param nthreads Number of processor threads for the filtering. If not
#'  specifed then the maximum number of logical cores are used.
#' @examples
#' Laurasmappings_filtered <- ZeroFilter(Laurasmappings, nthreads=4)
#'
#' @export

ZeroFilter <- function(dataset, non_zero_num = 4, nthreads = NULL) {
    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    genenumber <- nrow(dataset)
    timevector <- CircadianTools::MakeTimevector(dataset)
    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)
    filterdf <- foreach::foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        num_zeroes <- sum(as.list((dataset[i, ][-1])) == 0)  # Number of zeroes
        if (num_zeroes < (length(timevector) - non_zero_num)) {
            dataset[i, ]
        }
    }
    parallel::stopCluster(cl)
    rownames(filterdf) <- seq(1, nrow(filterdf))  #Rebuild the row names
    return(filterdf)
}


#' CombiFilter
#' @description Filters a transcriptomics dataset by using ZeroFilter,
#'  AnovaFilter and SizeFilter.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param non_zero_num The minimum number of non-zero activity readings for
#'  genes in the filtered dataset.
#' @param threshold Set the p-value threshold for the filtering.
#' @param cutoff Proportion of total number of genes (which show at least some
#'  activity) to be removed.
#' @param nthreads Number of processor threads for the filtering. If not
#'  specified then the maximum number of logical cores are used.
#' @param anovafilter Logical. If anovafiltering should be used. Defaults to
#'  TRUE.
#' @param zerofilter Logical. If zero filtering should be used. Defaults to
#'  TRUE.
#' @param sizefilter Logical. If size filtering should be used. Defaults to
#'  TRUE.
#' @examples
#' Laurasmappings_filtered <- CombiFilter(Laurasmappings, nthreads=4)
#'
#' @export
CombiFilter <- function(dataset, non_zero_num = 4, threshold = 0.05,
                        cutoff = 0.1, anovafilter = TRUE,
                        zerofilter = TRUE, sizefilter = TRUE, nthreads = NULL) {

    if (anovafilter == TRUE) {
        # Filter via ANOVA
        dataset <- CircadianTools::AnovaFilter(dataset = dataset,
                                               threshold = threshold,
                                               nthreads = nthreads)

        # Garbage collection to return RAM to OS before next function
        invisible(gc())
    }

    if (zerofilter == TRUE) {
        # Remove genes showing very little activity
        dataset <- CircadianTools::ZeroFilter(dataset = dataset,
                                              non_zero_num = non_zero_num,
                                              nthreads = nthreads)

        # Garbage collection to return RAM to OS before next function
        invisible(gc())
    }

    if (SizeFilter == TRUE) {
        # Filter by removing genes with the smallest range
        dataset <- CircadianTools::sizefilter(dataset = dataset,
                                              cutoff = cutoff,
                                              nthreads = nthreads)

        # Garbage collection to return RAM to OS before next function
        invisible(gc())
    }
    return(dataset)
}


#' TFilter:
#' @description Applies a filter where a t-test is carried out on
#'  gene activity levels between time points. The number of significant changes
#'  between time points is found. If there is a sufficient number of significant
#'  changes and close to as many positive changes as negative changes then the
#'  gene is included in the filtered dataset.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param maxdifference The maximum difference between the number of signifcant
#'  positive t statistic and the number of signifcant negative t statistic.
#' @param minchanges The minimum number of significant changes found via t-tests
#'  for the gene to be included in the filtered dataset.
#' @param nthreads Number of processor threads to be used for calculating the
#'  distance matrix. If not specified then the maximum number of logical cores
#'  are used.
#' @param psignificance The maximum p-value for which a result of a t-test is
#'  classed as significant.
#' @return Returns a filtered transcriptomics dataset
#' @examples
#' filterdf <- TFilter(Laurasmappings)
#'
#' @export

TFilter <- function(dataset, maxdifference = 1, minchanges = 2,
                    psignificance = 0.05, nthreads = NULL) {

    # Load the dopar binary operator from foreach package
    `%dopar%` <- foreach::`%dopar%`

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)



    filterdf <- foreach::foreach(i = 1:nrow(dataset),
                                 .combine = rbind) %dopar% {

        # TAnalysis returns two values as a vector. First values represents a
        # positive significant change between time points whilst second value
        # represents a negative significant change.
        ups.downs <- TAnalysis(row.no = i, dataset = dataset,
                               psignificance = psignificance)

      # Finds the difference between the number of positive and negative changes
        observed.difference <- abs(ups.downs[1] - ups.downs[2])
        # The total number of significant changes
        total.changes <- ups.downs[1] + ups.downs[2]

        if (total.changes >= minchanges) {

            if (observed.difference <= maxdifference) {
                # Add the gene as part of the filtered dataset
                data.frame(dataset[i, ])
            }
        }
    }

    parallel::stopCluster(cl)
    return(filterdf)
}

