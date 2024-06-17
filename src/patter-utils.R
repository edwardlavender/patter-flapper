#' @title {patter} helpers

# Convert .acoustics to glatos format
as_glatos <- function(.acoustics) {
  data.frame(detection_timestamp_utc = .acoustics$timestamp, 
             transmitter_codespace = "000",
             transmitter_id = as.character(.acoustics$individual_id), 
             receiver_sn = as.character(.acoustics$receiver_id)
  )
}