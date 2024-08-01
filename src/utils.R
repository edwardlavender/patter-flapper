date_to_POSIXct <- function(x) {
  as.POSIXct(paste(x, "00:00:00"), tz = "UTC")
}
