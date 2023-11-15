# Explicitly register relevant suggested {patter} dependencies with {renv}

if (!requireNamespace("sparklyr", quietly = TRUE)) {
  renv::install("sparklyr", prompt = FALSE)
}
if (!requireNamespace("spatialEco", quietly = TRUE)) {
  renv::install("spatialEco@2.0-1", prompt = FALSE)
}
if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
  renv::install("spatstat.geom", prompt = FALSE)
}
if (!requireNamespace("spatstat.explore", quietly = TRUE)) {
  renv::install("spatstat.explore", prompt = FALSE)
}
if (!requireNamespace("zoo", quietly = TRUE)) {
  renv::install("zoo", prompt = FALSE)
}
if (!requireNamespace("qs", quietly = TRUE)) {
  renv::install("qs", type = "source", dependencies = TRUE, prompt = FALSE)
}
if (!requireNamespace("zoo", quietly = TRUE)) {
  renv::install("zoo", prompt = FALSE)
}