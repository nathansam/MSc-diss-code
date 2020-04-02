#' DendogramPlot
#' @description Plots the dendogram for a cluster in a clustered dataset.
#' @param cluster.no The number which identifies the cluster. If a vector is
#'  provided then multiple clusters will be used.
#' @param cluster.dataset A transcriptomics dataset where the final column
#'  details the cluster the gene belongs to. First column should be gene names.
#'  All remaining columns should be expression levels.
#' @param method The clustering method to be used. Currently accepts 'agglom'
#'  and 'diana'.
#' @param metric The distance metric to be used to calculate the distances
#' between genes. See parallelDist::parDist for all accepted arguments. Also
#' allows the option of 'abs.correlation'.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#'  TRUE
#' @param save  Logical. If TRUE, saves plot to working directory. Defaults to
#'  FALSE.
#' @param path The directory to be used for saving plots to. Uses the working
#'  directory by default. Not used if save = FALSE
#' @examples
#' t.filter <-TFilter(Laurasmappings)
#' diana.df <- DianaClustering(t.filter, k= 95)
#' DendogramPlot(5, diana.df, method = 'diana')
#' @export
DendogramPlot <- function(cluster.no, cluster.dataset, method = "agglom",
                          metric = 'euclidean', print = TRUE, save = FALSE,
                          path = NULL ){

  if (save == TRUE) {
    if (is.null(path) == FALSE){
      if (dir.exists(path) == FALSE) {
  # If save==TRUE then create directory for saved plots if doesn't already exist
        dir.create(path)
      }
    }
  }

  dataset <- NULL
  for (i in cluster.no){
    sub <- subset(cluster.dataset, cluster == i) # Subset by cluster
    dataset <- rbind(dataset, sub)

  }


  dataset$cluster <- NULL
  # Calculate distance matrix
  distance <- CircadianTools::DistanceGen(dataset, metric=metric)
  attr(distance ,"Labels") <- 1:(nrow(dataset))


  ### Apply clustering and convert to dendogram object ###
  if (method == "agglom"){ # Agglomerative clustering
    hc <- as.dendrogram(hclust(distance))
  }
  if (method == "diana"){ # Divisive clustering using DIANA
    hc <- as.dendrogram(cluster::diana(distance))
  }


  ### Change shorthand names for methods/metrics to proper names ###
  if (method == "agglom"){
    method = "Agglomerative Clustering"
  }
  if (method == "diana"){
    method = "DIANA Clustering"
  }
  if (metric == "euclidean"){
    metric <- "Euclidean Distance"
  }
  if (metric == "abs.correlation"){
    metric <- "Absolute Pearson's Correlation"
  }

  ### Plotting ###
  # Create ggplot object
  p <- ggdendro::ggdendrogram(hc, rotate = FALSE)
  # Add title
  p <- p + ggplot2::ggtitle(paste("Cluster ", cluster.no," for ", method,
                                  " When  Using ", metric,
                                  " as the Distance Measure", sep=""))
  # Right-align title
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
  # Change text size
  p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))

  # Create filename
  filename <- paste(cluster.no, "_", deparse(substitute(cluster.dataset)),
                    ".png", sep = "")

  if (save == TRUE){
    ggplot2::ggsave(filename,
                    p, path = path, width = 10, height = 4.5, units = "in")
  }
  if (print == TRUE){
    return(p)
  }
}


#' DendogramDatasetPlot
#' @description Plots the dendogram for every cluster in a clustered dataset.
#' @param cluster.dataset A transcriptomics dataset where the final column
#'  details the cluster the gene belongs to. First column should be gene names.
#'  All remaining columns should be expression levels.
#' @param method The clustering method to be used. Currently accepts 'agglom'
#'  and 'diana'.
#' @param metric The distance metric to be used to calculate the distances
#' between genes. See parallelDist::parDist for all accepted arguments. Also
#' allows the option of 'abs.correlation'.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#'  TRUE
#' @param save  Logical. If TRUE, saves plot to working directory. Defaults to
#'  FALSE.
#' @param path The directory to be used for saving plots to. Uses the working
#'  directory by default. Not used if save = FALSE
#' @examples
#' t.filter <-TFilter(Laurasmappings)
#' diana.df <- DianaClustering(t.filter, k= 95)
#' DendogramDatasetPlot(diana.df, method = 'diana')
#' @export
DendogramDatasetPlot <- function(cluster.dataset, method = "agglom",
                          metric = 'euclidean', print = TRUE, save = FALSE,
                          path = NULL ){

  clusters <- unique(cluster.dataset$cluster)

  for (i in clusters){
    CircadianTools::DendogramPlot(cluster.no = i,
                                  cluster.dataset = cluster.dataset,
                                  method = method, metric = metric,
                                  print = print, save = save, path = path)
  }
}


#' DatasetDendogram
#' @param cluster.dataset A transcriptomics dataset where the final column
#'  details the cluster the gene belongs to. First column should be gene names.
#'  All remaining columns should be expression levels.
#' @param method The clustering method to be used. Currently accepts 'agglom'
#'  and 'diana'.
#' @param metric The distance metric to be used to calculate the distances
#' between genes. See parallelDist::parDist for all accepted arguments. Also
#' allows the option of 'abs.correlation'.
#' @param nthreads Number of processor threads to be used for calculating the
#'  distance matrix. If not specified then the maximum number of logical cores
#'  is used.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#'  TRUE
#' @param save  Logical. If TRUE, saves plot to working directory. Defaults to
#'  FALSE.
#' @param path The directory to be used for saving plots to. Uses the working
#'  directory by default. Not used if save = FALSE
#' @examples
#' t.filter <-TFilter(Laurasmappings)
#' diana.df <- DianaClustering(t.filter, k= 95)
#' DendogramPlot(diana.df, method = 'diana')
#' @export
DatasetDendogram <- function(cluster.dataset, method = "agglom",
                             metric = "euclidean", nthreads = NULL,
                             print = TRUE, save = FALSE, path = NULL){


  if (save == TRUE & is.null(path) == FALSE & dir.exists(path) == FALSE) {
  # Create directory for saved plots if doesn't already exist
        dir.create(path)
      }

  if (is.null(nthreads) == TRUE) {
    # Set the threads to maximum if none is specified
    nthreads <- parallel::detectCores()
  }

      # Load the dopar binary operator from foreach package
      `%dopar%` <- foreach::`%dopar%`
      cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
      doParallel::registerDoParallel(cl)

  unique.clusters <- unique(cluster.dataset$cluster)  # List of cluster labels

  # Mean time profile for each cluster
  mean.clusters <- foreach::foreach(i = unique.clusters,
                                    .combine = rbind) %dopar% {
                     ClusterTimeProfile(i, cluster.dataset = cluster.dataset,
                                        nthreads = 1)
                                    }
  parallel::stopCluster(cl)

  row.names(mean.clusters) <- unique.clusters

  # Generate distance matrix
  if (metric == "abs.correlation"){
    distance <- CircadianTools::AbsCorDist(mean.clusters)
  } else {
    distance <- parallelDist::parDist(mean.clusters, method = metric,
                                      threads = nthreads)
  }

  if (method == "agglom"){ # Agglomerative clustering
    hc <- as.dendrogram(hclust(distance))
  }
  if (method == "diana"){ # Divisive clustering using DIANA
    hc <- as.dendrogram(cluster::diana(distance))
  }

  ### Change shorthand names for methods/metrics to proper names ###
  if (method == "agglom"){
    method = "Agglomerative Clustering"
  }
  if (method == "diana"){
    method = "DIANA Clustering"
  }
  if (metric == "euclidean"){
    metric <- "Euclidean Distance"
  }
  if (metric == "abs.correlation"){
    metric <- "Absolute Pearson's Correlation"
  }

  ### Plotting ###
  p <- ggdendro::ggdendrogram(hc, rotate = FALSE)
  # Add title
  p <- p + ggplot2::ggtitle(paste( method, " When  Using ", metric,
                                  " as the Distance Measure", sep=""))
  # Right-align title
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
  # Change text size
  p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))

  # Create filename
  filename <- paste(deparse(substitute(cluster.dataset)),"_dendogram", ".png",
                    sep = "")

  if (save == TRUE){
    ggplot2::ggsave(filename,
                    p, path = path, width = 10, height = 4.5, units = "in")
  }
  if (print == TRUE){
    return(p)
  }
}





