#' Internal patter functions
if (!patter:::os_linux()) {
  utils.add::load_internal_functions("patter")
} else {
  nothing <- patter:::nothing
  msg     <- patter:::msg
  warn    <- patter:::warn
  abort   <- patter:::abort
  julia_check_exists <- patter:::julia_check_exists
}

#' glue
glue <- glue::glue
