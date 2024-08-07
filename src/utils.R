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

