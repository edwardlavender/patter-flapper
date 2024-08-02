rm(list = ls())
dv::src()

library(collapse)
library(patter)
library(lubridate)

#### Load data 
ud_grid   <- terra::rast(here_data("spatial", "ud-grid.tif"))
iteration <- qs::qread(here_data("input", "iteration", "coa.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))

#### Estimate coordinates (COAs)
terra:::readAll(ud_grid)
iteration_ls <- split(iteration, seq_row(iteration))[1:2]
cl_lapply(iteration_ls, function(sim) {
  # sim <- iteration[1, ]
  estimate_coord_coa(sim = sim, 
                     map = ud_grid,
                     datasets = list(detections = acoustics_by_unit[[sim$unit_id]],
                                     moorings = moorings))
})

#### Estimate UDs
estimate_uds(iteration)

