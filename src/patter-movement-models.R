#' @title Uniform movement models 

rlenuniform <- function(n, .mobility) {
  runif(n, min = 0, max = .mobility)
}

dlenuniform <- function(x, .mobility) {
  dunif(x, min = 0, max = .mobility)
}