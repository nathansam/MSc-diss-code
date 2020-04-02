#' MDSPlot
#' @description Applies multidimensional scaling to a clustered transcriptomics
#'  dataset to reduce the clusters to two dimensions and then plots the clusters. 
#' @param cluster.dataset A transcriptomics dataset where the final column
#'  details the cluster the gene belongs to. First column should be gene names.
#'  All remaining columns should be expression levels.
#' @param title The title to be used in the plot
#' @param nthreads Number of processor threads for the process. If not
#'  specified then the maximum number of logical cores are used.
#' @param metric The distance metric to be used to calculate the distances
#'  between genes. See parallelDist::parDist for all accepted arguments. Also 
#'  accepts "abs.correlation" for absolute Pearson's correlation
#' @param save Logical. If TRUE, saves plots. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer.
#'  Defaults to TRUE
#' @examples 
#' a.filter <- AnovaFilter(Laurasmappings)
#' cluster.df <- PamClustering(a.filter, k = 95, scale = TRUE)
#' MDSPlot(cluster.df)
#' @return 
#' 
MDSPlot <- function(cluster.dataset, title = "Multidimensional Scaling Plot",
                    nthreads = NULL, metric = "euclidean", save = FALSE,
                    print = TRUE){
  
  
  # Get the center for each cluster
  center <- CircadianTools::ClusterCenterGenerator(cluster.dataset = 
                                                     cluster.dataset,
                                                   nthreads = nthreads)
  
  # Calculate distance matrix 
  if (metric == "abs.correlation"){
    distance <- CircadianTools::AbsCorDist(center, nthreads = nthreads)
  } else {
    distance <- parallelDist::parDist(center, method = metric)
  }
  
  
  # Perform Multidimensional scaling
  fit <- cmdscale(distance, k = 2, eig = TRUE)
  x <- fit$points[,1] # First dimension
  y <- fit$points[,2] # Second dimension
  

  # Colours to be used for the plot in order to make points more readable
  colours.vector <- c("#412d6b", "#008dd5", "#ffa630", "#ba1200", "#840032")
  
  # Create dataframe of x,y coords and which cluster the row is for
  df <- data.frame(x, y, cluster = row.names(center))
  
  # Create ggplot opject
  p <- ggplot2::ggplot(ggplot2::aes(x = x, y = y, color = as.numeric(cluster)),
                       data = df)
  # Plots clusters as text
  p <- p + ggplot2::geom_text(ggplot2::aes(label = cluster),
                              show.legend = FALSE)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("Dimension 1")
  p <- p + ggplot2::ylab("Dimension 2")
  p <- p + ggplot2::ggtitle(title) # Either user provided or the default
  p <- p +ggplot2::scale_colour_gradientn(colors = colours.vector ) # Txt colour
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
  p <- p + ggplot2::theme(text = ggplot2::element_text(size = 12))

  
  if (save == TRUE){
    filename <- paste (deparse(substitute(cluster.dataset)), "_MDS_plot.png",
                       sep = "")
    ggplot2::ggsave(filename, plot = p, width = 10, height = 4.5,
                    units = "in")  # Save the plot
  }
  
  if (print == TRUE){
    return(p)
  }
}