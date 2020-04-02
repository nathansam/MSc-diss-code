#' ContigGen
#' @description Finds all unique contig IDs in a transcriptomics dataset
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#' All other columns should be expression levels.
#' @return A character vector of contig IDs
#' @examples
#' contigs <- ContigGen(Laurasmappings)
#'
#' @export
ContigGen <- function (dataset){
  contigs <- dataset$sample # Get sample IDs
  # Split strings by underscores into columns
  contigs <- read.table(text = contigs, sep = "_")

  # Join the characters before the first underscore and after the last
  contigs <- paste (contigs[,1], contigs[,3], sep = "_")

  contigs <-unique(contigs)
  return(contigs)
}

#' FastaSub
#' @description  Creates a fasta file from only certain sequences in another
#'  fasta file
#' @param gene.names A character vector of gene names
#' @param fasta.file A fasta file read in using \code{seqinr}
#' @param save Logical. If TRUE then the newly created fasta file is saved to
#'  the working directory
#' @param filename The filename to use for a saved fasta file
#' @return A fasta file
#' @examples
#' main.fasta <- seqinr::read.fasta("~/JoesTranscriptomeMin300bp2.fasta")
#' a.filter <- AnovaFilter(Laurasmappings)
#' fasta.sub <- FastaSub(a.filter$sample, main.fasta, filename="a_filter")
#' @export
 FastaSub <- function(gene.names, fasta.file, save=TRUE, filename=NULL){

   if (is.null(filename)==TRUE){
     filename <- deparse(substitute(fasta.file))
     filename <- paste(filename, "_subset", sep = "")
   }

   filename <- paste(filename, ".fasta", sep="")

   contigs <- unique(gene.names)
   fasta.sub <- fasta.file[names(fasta.file) %in% contigs]
   seqinr::write.fasta(fasta.sub,names= names(fasta.sub), file.out = filename)

   return(fasta.sub)
 }
