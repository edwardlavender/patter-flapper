#' Internal patter functions
if (!patter:::os_linux()) {
  utils.add::load_internal_functions("patter")
} else {
  
  check_timeline   <- patter:::check_timeline
  check_named_list <- patter:::check_named_list
  check_names      <- patter:::check_names

  julia_check_exists <- patter:::julia_check_exists
  
  msg     <- patter:::msg
  warn    <- patter:::warn
  abort   <- patter:::abort
  
  nothing <- patter:::nothing
  
}

#' glue
glue <- glue::glue
