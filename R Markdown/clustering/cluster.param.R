library(CircadianTools)


a.filter <- read.csv(
  "~/MEGA/Uni/Masters/Diss/Stats/R Markdown/filtering/a_filter.csv",
  stringsAsFactors=FALSE)
a.filter <-GeneScale(a.filter)
t.filter <- read.csv(
  "~/MEGA/Uni/Masters/Diss/Stats/R Markdown/filtering/t_filter.csv",
  stringsAsFactors=FALSE)
t.filter <-GeneScale(t.filter)

set.seed(123)

k.params <- c(seq(2,5), seq(10, 300, 5))


nthreads <- parallel::detectCores()
  
  ClusterParamSelection(a.filter, k = k.params, metric="euclidean",
                        nthreads = nthreads, path = "anova_euclid_validation")


ClusterParamSelection(a.filter, k = k.params, metric="abs.correlation",
                      nthreads = nthreads, path = "anova_abscor_validation")


ClusterParamSelection(t.filter, k = k.params, metric="euclidean",
                      nthreads = nthreads, path = "ttest_euclid_validation")

ClusterParamSelection(t.filter, k = k.params, metric="abs.correlation",
                      nthreads = nthreads, path = "ttest_abscor_validation")