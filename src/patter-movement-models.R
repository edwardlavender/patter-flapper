#' @title Uniform movement models 

rlenuniform <- function(.n, .mobility) {
  runif(.n, min = 0, max = .mobility)
}

dlenuniform <- function(.x, .mobility) {
  dunif(.x, min = 0, max = .mobility)
}

rkickuniform <- function(.xy0,
                         .rlen = rlenuniform,
                         .rang = rangrw,
                         .obs, .t, .dlist) {
  # Simulate step lengths and turning angles
  n    <- nrow(.xy0)
  rlen <- .rlen(n, .mobility = .obs$mobility[.t])
  rang <- .rang(n)
  # Step into new locations
  cstep(.xy0 = .xy0,
        .len = rlen,
        .ang = rang,
        .lonlat = .dlist$pars$lonlat)
}

dkickuniform <- function(.xy0, .xy1, .obs, .t, .dlist) {
  .xy0 <- as.matrix(.xy0)
  .xy1 <- as.matrix(.xy1)
  # Calculate step length between selected location and all previous locations
  rlen <- clen(.xy0 = .xy0, .xy1 = .xy1, .lonlat = .dlist$pars$lonlat)
  # Calculate densities 
  dlenuniform(rlen, .mobility = .obs$mobility[.t])
}

#' @title Truncated gamma movement models 

rkickgamma <- function(.xy0,
                         .rlen = rlen,
                         .rang = rangrw,
                         .obs, .t, .dlist) {
  # Simulate step lengths and turning angles
  n    <- nrow(.xy0)
  rlen <- .rlen(n, .shape = .obs$shape[.t], .scale = .obs$scale[.t], .mobility = .obs$mobility[.t])
  rang <- .rang(n)
  # Step into new locations
  cstep(.xy0 = .xy0,
        .len = rlen,
        .ang = rang,
        .lonlat = .dlist$pars$lonlat)
}

dkickgamma <- function(.xy0, .xy1, .obs, .t, .dlist) {
  .xy0 <- as.matrix(.xy0)
  .xy1 <- as.matrix(.xy1)
  # Calculate step length between selected location and all previous locations
  rlen <- clen(.xy0 = .xy0, .xy1 = .xy1, .lonlat = .dlist$pars$lonlat)
  # Calculate densities 
  dtruncgamma(rlen, .shape = .obs$shape[.t], .scale = .obs$scale[.t], .mobility = .obs$mobility[.t])
}
