#' @title {Tools4ETS} imports & extensions

# Differences
difference <- function(x2, x1, f = NULL, ...) {
  if (class(x2)[1] %in% c("numeric", "integer")) {
    d <- x2 - x1
  }
  else if (class(x2)[1] %in% c("POSIXct", "POSIXlt", "Date")) {
    d <- difftime(x2, x1, ...)
  }
  if (!is.null(f)) {
    d <- as.numeric(d)
  }
  d
}

# Serial differences 
serial_difference <- function(x, na.rm = FALSE, ...) {
  dur <- difference(dplyr::lead(x), x, ...)
  if (na.rm) {
    posNA <- which(is.na(dur))
    dur <- dur[-c(posNA)]
  }
  dur
}

# Month/year category ranges
mmyyrng <- function(mmyy) {
  stopifnot(length(mmyy) == 1L)
  start  <- lubridate::dmy(paste0("01-", mmyy))
  start  <- paste(start, "00:00:00")
  start <- as.POSIXct(start, tz = "UTC")
  end   <- start + lubridate::period(1, "month") - lubridate::period("2 mins")
  c(start, end)
}

# Examples 
if (FALSE) {
  lapply(paste0("0", 1:9, "-2016"), mmyyrng)
  lapply(paste0(10:12, "-2016"), mmyyrng)
}