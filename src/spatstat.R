#' simplify.owin.trial()

# Smoothing the window improves UD estimation speed
if (!patter:::os_linux()) {

# We visualise different window smoothnesses
# And test how long UD fitting takes with a sample of data
# This code is based on https://github.com/spatstat/spatstat.explore/issues/4
simplify.owin.trial <- function(sea, win, dmin = 0) {
  
  #### Define example points
  xy <- 
    terra::spatSample(sea, size = 200L) |> 
    terra::crds() |> 
    as.data.frame()

  #### Simplify window by dmin
  win_s <- win
  if (dmin > 0) {
    cat("Simplifying window...\n")
    tic()
    win_s <- spatstat.geom::simplify.owin(win, dmin)
    toc()
  }
  
  #### Visualise window 
  cat("Visualising window...\n")
  tic()
  png(here_fig("spatstat", "windows", paste0(dmin, "-window.png")), 
      height = 10, width = 10, units = "in", res = 600)
  plot(win)
  plot(win_s, add = TRUE, border = "green")
  dev.off()
  toc()
  
  #### Define ppp
  cat("Defining ppp using win_s...\n")
  X <- spatstat.geom::ppp(x = xy$x, y = xy$y, window = win_s)
  
  #### Estimate UD (500 pixels)
  cat("Estimating UD...")
  spatstat.geom::spatstat.options(npixel = 500)
  tic()
  D <- spatstat.explore::density.ppp(X, at = "pixels")
  toc()
  
  #### Plot UD
  cat("Plotting UD...\n")
  tic()
  png(here_fig("spatstat", "windows", paste0(dmin, "-ud.png")), 
      height = 10, width = 10, units = "in", res = 600)
  plot(D)
  dev.off()
  toc()

  win_s
  
}

}
