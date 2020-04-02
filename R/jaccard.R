#' Jaccard
#'
#' @description Stolen from the clusteval package as there is currently a bug
#'  with the package on CRAN where this function is not exported from the
#'  namespae of clusteval. Computes the Jaccard similarity coefficient of two
#'  clusterings of the same data set under the assumption that the two
#'  clusterings are independent.
#'
#'  For two clusterings of the same data set, this function calculates the
#'  Jaccard similarity coefficient of the clusterings from the comemberships of
#'  the observations. Basically, the comembership is defined as the pairs of
#'  observations that are clustered together.
#'

#' @details
#'  The Jaccard similarity coefficient is defined as:
#'  \deqn{J = \frac{n_{11}}{n_{11} + n_{10} + n_{01}}}.
#'
#'  In the special case that the Jaccard coefficient results in \eqn{0/0},
#'  we define \eqn{J = 0}. For instance, this case can occur when both clusterings
#'  consist of all singleton clusters.
#'
#'  To compute the contingency table, we use the \code{\link{comembership_table}}
#'  function.
#'
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Jaccard coefficient for the two sets of cluster labels (See
#' Details.)
#' @examples
#
#'
#' @export
Jaccard <- function(labels1, labels2) {
    com_table <- clusteval::comembership_table(labels1, labels2)
    jaccard_out <- with(com_table, n_11/(n_11 + n_10 + n_01))

    # In the case where 'labels1' and 'labels2' contain all singletons, the Jaccard coefficient results in
    # the expression 0 / 0, which yields a NaN value in R.  We define such cases as 0.
    if (is.nan(jaccard_out)) {
        warning("The two clusterings contain all singletons -- returning 0.")
        jaccard_out <- 0
    }
    jaccard_out
}
