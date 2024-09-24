#' Platform wrappers

on_unix <- function() {
  .Platform$OS.type == "unix"
}

#' Cat helpers

cat_line <- function() {
  cat("\n\n\n---------------------------------------------------------------\n")
}

cat_iter <- function(t) {
  msg("\n On row {t}...", .envir = environment())
}

cat_init <- function(t) {
  cat_line()
  cat_iter(t)
}

#' Create directories 
dirs.create <- function(paths) {
  sapply(paths, function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
  }) |> 
    suppressWarnings()
  invisible(NULL)
}

#' Sort bathymetry datasets

sort_bathysets <- function(x) {
  x |> dirname() |> basename() |> unique() |> gtools::mixedsort()
}

fill <- function(x, name) {
  c(name, rep("", length(x) - 1))
}