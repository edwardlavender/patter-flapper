#' Internal patter functions
if (!patter:::os_linux()) {
  utils.add::load_internal_functions("patter")
} else {
  nothing <- patter:::nothing
}


#' glue
glue <- glue::glue
