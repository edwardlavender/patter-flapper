progress <- function(t, nt = max(t)) {
  msg("On row {t} / {nt} ({round(t/nt * 100)} %)...", .envir = environment())
}

dirs.create <- function(paths) {
  sapply(paths, function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
  }) |> 
    suppressWarnings()
  invisible(NULL)
}

