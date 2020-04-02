#' MultiCosinorTest
#' @description Fits a cosinor model and carries out ANOVA using raw
#'  coefficients. Then fits a cosinor model with additonal sine and cosine terms
#'  with a different period. ANOVA tests are carried out on the more complex
#'  model as well as directly comparing the two models.
#'
#'
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'  All other columns should be expression levels.
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param period.1 The period of the simple Cosinor Model.
#' @param period.2 The period of the additional cosinor terms in the more
#'  complex model.
#' @param save Logical. If TRUE then the ANOVA test results will be saved to
#'  three csv files.
#' @param path The directory to be used for saving ANOVA results to. Only used
#'  if save == TRUE Uses theprovided genename appended with '_cosinor_tests' if
#'  this argument is not specified.
#' @examples
#' MultiCosinorTest("comp99801_c1_seq1", Laurasmappings)
#' @export
MultiCosinorTest <- function(genename, dataset, timelag = 6, period.1 = 24,
                          period.2 = 12.4, save = TRUE, path = NULL ){

  if (save == TRUE){
    if (is.null(path) == TRUE){
      path = paste(genename,"_cosinor_tests", sep="" )
    }

    if (dir.exists(path) == FALSE) {
      # Create directory for saved plots if it doesn't already exist
      dir.create(path)
    }
  }

  activity <- subset(dataset, dataset[1] == genename)
  timevector <- CircadianTools::MakeTimevector(activity)  # Makes timevector

  activity <- activity[-1]  # Remove gene name from subset
  activity <- t(activity)
  data <- data.frame(timevector - timelag, activity)


  names(data) <- c("timevector", "activity")

  data$rrr.1 <- cos((2 * pi * timevector) / period.1)
  data$sss.1 <- sin((2 * pi * timevector) / period.1)

  data$rrr.2 <- cos((2 * pi * timevector) / period.2)
  data$sss.2 <- sin((2 * pi * timevector) / period.2)

  simple <- lm(activity ~ rrr.1 + sss.1, data = data)
  complex <- lm (activity ~ rrr.1 + sss.1 + rrr.2 + sss.2, data = data)

  cat(crayon::red("The anova table for the simple model is given by: \n"))
  simple.anova <- anova(simple)
  print(simple.anova)

  cat(crayon::red("The anova table for the more complex model is given by: \n"))
  complex.anova <- anova(complex)
  print(complex.anova)

  cat(crayon::red("Carrying out ANOVA on both models produces the following table : \n"))
  comparison.anova <- anova (simple, complex)
  print(comparison.anova)

  if (save == TRUE){
    write.csv(simple.anova, file = paste(path,"/simple.csv", sep=""))
    write.csv(complex.anova, file = paste(path,"/complex.csv", sep =""))
    write.csv(comparison.anova, file = paste(path,"/comparison.csv", sep=""))
    cat(crayon::red("The ANOVA results have been saved. \n"))
  }

}
