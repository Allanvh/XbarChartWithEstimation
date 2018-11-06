#' prepare_package
#'
#'
#' @import cubature, numDeriv
#'
#'@export


prepare_package <- function () {
install.packages("cubature")
library(cubature)

install.packages("numDeriv")
library(numDeriv)

}
