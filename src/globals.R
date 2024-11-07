#' Internal patter functions
if (!patter:::os_linux()) {
  utils.add::load_internal_functions("patter")
} else {
  
}


#' glue
glue <- glue::glue
