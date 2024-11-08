#' Internal patter functions
if (!patter:::os_linux()) {
  utils.add::load_internal_functions("patter")
} else {
  nothing <- patter:::nothing
  msg     <- patter:::msg
  warn    <- patter:::warn
  abort   <- patter:::abort
}


#' glue
glue <- glue::glue
