
estimate_coord_coa <- function(sim, map, datasets) {
  
  # Estimate COAs
  t1 <- Sys.time()
  coord <- coa(.map = map,
               .acoustics = datasets$detections,
               .moorings = datasets$moorings,
               .delta_t = sim$delta_t,
               .plot_weights = FALSE)
  t2 <- Sys.time()
  
  # Record convergence & timing
  time <- sim[, time := secs(t2, t1)]
  
  # Save outputs
  qs::qsave(time, file.path(sim$folder, "time.qs"))
  qs::qsave(coord, file.path(sim$folder, "coord.qs"))
  nothing()
  
}

