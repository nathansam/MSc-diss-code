#' DistanceGen
#' @description Generates a distance matrix from a a transcriptomics dataset.
#' @param dataset A transcriptomics dataset. Preferably filtered first. First
#' columns should be gene names. All other columns should be expression levels.
#' @param metric The distance metric to be used to calculate the distances
#'  between genes. See parallelDist::parDist for all accepted arguments. Also
#'  allows the option of 'abs.correlation'. Not used if a distance matrix is
#'  provided.
#' @param nthreads The number of threads to be used for parallel computations.
#'  If NULL then the maximum number of threads available will be used.
#' @examples
#' a.filter <- AnovaFilter(Laurasmappings)
#' distance <- DistanceGen(a.filter, metric='abs.correlation')
#'
#' @export
DistanceGen <- function(dataset, metric = "euclidean", nthreads = NULL) {
    # Calculate the medians at each timepoint
    dataset <- CircadianTools::MedList(dataset, nthreads = nthreads)
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }
    if (metric == "abs.correlation") {
        distance <- AbsCorDist(dataset)
    } else{
        #Calculate the distance matrix
        distance <- parallelDist::parDist(dataset, method = metric,
                                          threads = nthreads)
    }
    return(distance)
}
