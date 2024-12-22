###########################
###########################
#### run-patter.R

#### Aims
# 1) Run patter algorithms

#### Prerequisites
# 1) Previous scripts


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
Sys.setenv(JULIA_SESSION = "FALSE")
dv::src()

#### Load data 
if (!patter:::os_linux() | Sys.getenv("JULIA_SESSION") == "FALSE") {
  map   <- terra::rast(here_data("spatial", "bathy.tif"))
  coast <- qreadvect(here_data("spatial", "coast.qs"))
}
pars              <- qs::qread(here_data("input", "pars.qs"))
iteration         <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))
col_sensitivity   <- readRDS(here_data("graphics", "col_sensitivity.rds"))


###########################
###########################
#### Run algorithm 

#### Julia Set up
if (patter:::os_linux() & Sys.getenv("JULIA_SESSION") == "TRUE") {
  stopifnot(!any(c("terra", "sf") %in% sort(loadedNamespaces())))
}
if (Sys.getenv("JULIA_SESSION") == "TRUE") {
  julia_connect()
  set_seed()
  set_map(here_data("spatial", "bathy.tif"))
  set_model_move_components()
}

#### Clean up
if (FALSE) {
  unlink(iteration$folder_coord, recursive = TRUE)
  unlink(iteration$folder_ud, recursive = TRUE)
  list.files(iteration$folder_coord, recursive = TRUE)
  list.files(iteration$folder_ud, recursive = TRUE)
  dirs.create(iteration$folder_coord)
  dirs.create(iteration$folder_ud)
  dirs.create(file.path(iteration$folder_ud, "spatstat", "h"))
}

#### Prepare iterations/datasets
iteration[, folder_xinit := file.path("data", "input", "xinit", "1.5", individual_id, month_id)]
iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings           = moorings,
                 archival_by_unit   = archival_by_unit, 
                 behaviour_by_unit  = behaviour_by_unit)

#### Review inputs
# Check iterations
# * Each unit_id is a individual/month combination
# * For each unit_id, we have AC/DC/ACDC algorithm implementations
# * For each algorithm implementation, we have different parameterisations (sensitivity)
head(iteration)
# Check datasets 
# * unit_ids were defined as _all possible_ individual-month combinations
# * (see process-data-mefs.R)
# -> We have 154 possible unit ids (individual-month) combinations
ids  <- unique(iteration$individual_id)
mons <- mmyy(seq(as.Date("2016-04-01"), as.Date("2017-05-31"), by = "months"))
nrow(CJ(ids, mons))
# * We have 48 'realised' unit_ids 
# -> These are the individual/months for which we have observations
sort(unique(iteration$unit_id))
length(unique(iteration$unit_id))
# * We have 154 elements in each dataset (one for each _possible_ unit_id)
length(datasets$detections_by_unit)
length(datasets$archival_by_unit)  
length(datasets$behaviour_by_unit)
# * For all AC/ACDC unit_id, we need detections dataset
sapply(unique(iteration$unit_id[iteration$dataset %in% c("ac", "acdc")]), function(id) {
  !is.null(datasets$detections_by_unit[[id]])
}) |> all() |> stopifnot()
# * For all DC/ACDC unit_id, we need archival dataset
sapply(unique(iteration$unit_id[iteration$dataset %in% c("dc", "acdc")]), function(id) {
  !is.null(datasets$archival_by_unit[[id]])
}) |> all() |> stopifnot()
# * For all AC/DC/ACDC unit_id, we need behaviour dataset
sapply(unique(iteration$unit_id), function(id) {
  !is.null(datasets$behaviour_by_unit[[id]])
}) |> all() |> stopifnot()
# Check folder_xinit
all(file.exists(file.path(iteration$folder_xinit, "xinit-fwd.qs"))) |> stopifnot()
all(file.exists(file.path(iteration$folder_xinit, "xinit-bwd.qs"))) |> stopifnot()
# Check number of rows
table(iteration$dataset)
table(iteration$mobility)
table(iteration$dataset, iteration$mobility)

#### Estimate coordinates
if (FALSE) {
  
  #### Set batch/rows
  # Set batch & rows
  # * batch 1: 528 rows
  # * batch 2: 144 rows 
  # * batch 3: 144 rows
  cmdargs <- commandArgs(trailingOnly = TRUE)
  # batch <- pars$pmovement$mobility[1]
  batch <- as.numeric(cmdargs[1])
  rows  <- eval(parse(text = cmdargs[2]))
  # rows <- 1L
  # Select batch & rows
  iteration <- iteration[mobility == batch, ]
  iteration <- iteration[rows, ]
  # Export vmap
  vmap      <- here_data("spatial", glue("vmap-{batch}.tif"))
  set_vmap(.vmap = vmap)
  rm(vmap, cmdargs)
  
  #### Clean memory 
  if (TRUE) {
    # pryr::mem_used() # 314 MB
    objs      <- ls()
    objs      <- objs[sapply(objs, function(x) !is.function(get(x)))]
    objs_keep <- c("iteration", "datasets", 
                   "batch", "rows", 
                   "logpdf_step.ModelMoveFlapper", "logpdf_step.ModelMoveFlapperCRW", 
                   "ModelMoveFlapperCRW", "moorings", "rho", 
                   "simulate_step.ModelMoveFlapper", "simulate_step.ModelMoveFlapperCRW", 
                   "state_flapper", "states_init.StateXYD")
    rm(list = objs[!(objs %in% objs_keep)])
    # pryr::mem_used() # 314 MB
  }
  
  #### Estimate coordinates
  gc()
  nrow(iteration)
  lapply_estimate_coord_patter(iteration  = iteration,
                               datasets   = datasets, 
                               trial      = FALSE, 
                               log.folder = here_data("output", "log", "real"), 
                               log.txt    = glue("log-{batch}-{min(rows)}.txt"))
  
}

#### Patter error check 
# (Code copied from simulate-algorithms.R)
if (FALSE) {
  
  # Check for errors on forward filter: OK 
  sapply(split(iteration, seq_row(iteration)), function(d) {
    qs::qread(file.path(d$folder_coord, "data-fwd.qs"))$error
  }) |> unlist() |> unique()
  
  # Check for errors on backward filter: OK
  sapply(split(iteration, seq_row(iteration)), function(d) {
    file <- file.path(d$folder_coord, "data-bwd.qs")
    if (file.exists(file)) {
      qs::qread(file)$error
    }
  }) |> unlist() |> unique()
  
  # Check for errors on smoother: OK
  sapply(split(iteration, seq_row(iteration)), function(d) {
    file <- file.path(d$folder_coord, "data-smo.qs")
    if (file.exists(file)) {
      qs::qread(file)$error
    }
  }) |> unlist() |> unique()
  
  # Check file sizes (MB) for reference
  # > Smoothed particles for one month are ~930 MB
  # > (2000 smoothing particles)
  sapply(split(iteration, seq_row(iteration)), function(d) {
    file <- file.path(d$folder_coord, "coord-smo.qs")
    if (file.exists(file)) {
      file.size(file) / 1e6
    }
  }) |> unlist() |> utils.add::basic_stats()
  
  # min   mean median    max    sd   IQR   MAD
  # 742.87 920.63 930.34 990.37 43.34 48.53 35.83
  
}

#### Estimate UDs
if (TRUE && (!patter:::os_linux() | Sys.getenv("JULIA_SESSION") == "FALSE")) {
  
  # Quick convergence check
  check <- copy(iteration)
  check[, convergence := file.exists(file_coord)][, .(index, convergence)] |> 
    as.data.frame()
  table(check$convergence)
  
  # Examine selected coord
  # lapply_qplot_coord(iteration, 
  #                    "coord-smo.qs",
  #                    extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])
  
  #### Estimate UDs
  # Time trial 
  # * 26 s: siam-linux20
  lapply_estimate_ud_spatstat(iteration     = iteration[1, ], 
                              extract_coord = function(s) s$states,
                              cl            = NULL, 
                              plot          = FALSE)
  # Implementation 
  # * (~43 min)
  lapply_estimate_ud_spatstat(iteration     = iteration, 
                              extract_coord = function(s) s$states,
                              cl            = 10L, 
                              plot          = FALSE)
  # (optional) Examine selected UDs
  lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
  
}


###########################
###########################
#### Generate outputs


###########################
#### Tidy iteration (copied from simulate-algorithms.R)

if (TRUE) {
  
  # Define grouping variables with tidy labels
  iteration[, unit_id := factor(unit_id)]
  iteration[, dataset := factor(dataset, levels = c("ac", "dc", "acdc"), labels = c("AC", "DC", "ACDC"))]
  # Revise 'sensitivity' coding
  iteration[sensitivity == "movement" & mobility == pars$pmovement$mobility[2], sensitivity := "movement-under"]
  iteration[sensitivity == "movement" & mobility == pars$pmovement$mobility[3], sensitivity := "movement-over"]
  iteration[sensitivity == "ac" & receiver_alpha == pars$pdetection$receiver_alpha[2], sensitivity := "ac-under"]
  iteration[sensitivity == "ac" & receiver_alpha == pars$pdetection$receiver_alpha[3], sensitivity := "ac-over"]
  iteration[sensitivity == "dc" & depth_sigma == pars$pdepth$depth_sigma[2], sensitivity := "dc-under"]
  iteration[sensitivity == "dc" & depth_sigma == pars$pdepth$depth_sigma[3], sensitivity := "dc-over"]
  iteration[, sensitivity := factor(sensitivity, 
                                    levels = c("best", 
                                               "movement-under", "movement-over", 
                                               "ac-under", "ac-over",
                                               "dc-under", "dc-over"), 
                                    labels = c("Best", 
                                               "Move (-)", "Move (+)", 
                                               "AC(-)", "AC(+)", 
                                               "DC(-)", "DC(+)"))] 
  qs::qsave(iteration, here_data("input", "iteration", "patter-tidy.qs"))
  
}


###########################
#### Convergence 

if (TRUE) {
  
  # Determine convergence
  iteration[, convergence := sapply(split(iteration, seq_row(iteration)), function(d) {
    file.exists(file.path(d$folder_coord, "coord-smo.qs"))
  })]
  
  # Summarise convergence
  iteration |> 
    group_by(dataset, sensitivity) |> 
    summarise(convergence = length(which(convergence)) / n()) |>
    ungroup() |>
    group_by(dataset) |> 
    mutate(range = max(convergence) - min(convergence)) |>
    ungroup() |>
    arrange(dataset, convergence) |>
    as.data.table()
  
  # Plot convergence proportions 
  png(here_fig("analysis", "convergence.png"), 
      height = 3, width = 9, units = "in", res = 600)
  iteration |> 
    group_by(dataset, sensitivity) |> 
    summarise(convergence = length(which(convergence)) / n()) |>
    ungroup() |>
    # tidyr::complete(dataset, sensitivity) |>
    as_tibble() |> 
    ggplot(aes(dataset, convergence, fill = sensitivity)) +
    geom_bar(stat = "identity", 
             position = position_dodge(preserve = "single"),
             colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = col_sensitivity) + 
    # scale_x_discrete(drop = FALSE) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
    labs(
      x = "Algorithm",
      y = "Pr(convergence)",
      fill = "Parameterisation"
    ) + 
    theme_bw() + 
    theme(axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          panel.grid.major.x = element_line(color = "gray80", linewidth = 0.2),
          panel.grid.minor.x = element_line(color = "gray80", linewidth = 0.2),
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank())
  dev.off()
  
}


###########################
#### Diagnostics
# (copied from simulate-algorithms.R)

if (TRUE) {
  
  #### Extract standard diagnostics (ESS) (~2 min)
  iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
  diagnostics <- 
    cl_lapply(iteration$file_coord, .fun = function(file_coord) {
      if (file.exists(file_coord)) {
        diag <- qs::qread(file_coord)$diagnostics
        diag[, file_coord := file_coord]
        diag[, .(file_coord, ess, maxlp)]
      }
    }, .cl = 10L, .combine = rbindlist)
  # Join iteration and diagnostics
  diagnostics <- full_join(iteration, diagnostics, by = "file_coord")
  
  #### Check smoothing success
  diagnostics |> 
    filter(convergence) |> 
    group_by(file_coord) |> 
    summarise(nan_perc = length(which(is.na(ess))) / n() * 100) |> 
    arrange(desc(nan_perc)) |> 
    as.data.table()
  
  #### Examine ESS
  # Visualise histogram of ESS for 'best' simulations 
  diagnostics |> 
    filter(convergence & sensitivity == "Best" & iter == 1L) |> 
    ggplot() + 
    geom_histogram(aes(ess), bins = 500) + 
    scale_y_continuous(expand = c(0, 0)) + 
    xlab("ESS") + ylab("Frequency") + 
    facet_wrap(~dataset) +
    theme_bw() + 
    theme(panel.grid = element_blank()) 
  # Compute summary statistics
  diagnostics |> 
    filter(convergence & sensitivity == "Best" & iter == 1L) |> 
    # group_by(dataset) |>
    summarise(utils.add::basic_stats(ess, na.rm = TRUE)) |>
    as.data.table()
  
}


###########################
#### Mapping (~1 min)

if (TRUE) {
  
  unique(iteration$dataset)
  unique(iteration$sensitivity)
  length(unique(iteration$individual_id))
  length(unique(iteration$month_id))
  length(unique(paste(iteration$sensitivity, iteration$dataset)))
  cl_lapply(unique(iteration$sensitivity), .cl = 9L, .fun = function(sens){
    lapply(c("AC", "DC", "ACDC"), function(alg) {
      
      # Define mapfiles
      mapfiles <-
        iteration |> 
        filter(sensitivity == sens & dataset == alg) |>
        mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif"), 
               individual_id = factor(individual_id, levels = sort(unique(individual_id)))) |>
        dplyr::select(row = individual_id, column = month_id, mapfile) |>
        as.data.table() 
      
      # Make map
      if (nrow(mapfiles)) {
        png_args <- 
          list(filename = here_fig("analysis", glue("map-patter-{alg}-{sens}.png")), 
               height = 10, width = 9, units = "in", res = 800)
        ggplot_maps(mapdt = mapfiles, png_args = png_args)
      }
      NULL
    })
  })
  
}


###########################
#### Geographic uncertainty

if (TRUE) {
  
  #### (1) Simplify coast for speed
  coast_s   <- terra::simplifyGeom(coast, tolerance = 100)
  coast_s_w <- terra::wrap(coast_s)
  
  #### (2) Test method
  # Sample points 
  pxy <- terra::spatSample(map, size = 100, xy = TRUE, na.rm = TRUE)
  # Compute MCP + cut out coast
  mcp <- terra::convHull(terra::vect(cbind(pxy$x, pxy$y), crs = terra::crs(map)))
  mcp <- terra::erase(mcp, coast_s)
  # Visual check
  terra::plot(map)
  points(pxy$x, pxy$y)
  terra::lines(mcp)
  # Compute area spanned by MCP
  terra::expanse(mcp, unit = "km")
  
  #### (3) Implementation
  # TO DO 
  # Rerun
  # Initial runtime: 
  # * 1 cl  : 42 min; 
  # * 2 cl  : 21 min
  # * 10 cl : 37 min
  # Limited benefits of cluster here due to memory requirements
  # (Updated run expected to take much longer due to larger files)
  tic()
  particle_mcps <- 
    cl_lapply(split(iteration, iteration$index), .cl = 2L, .fun = function(d) {
      # d <- iteration[1, ]
      print(d$index)
      if (file.exists(d$file_coord)) {
        # Read particles 
        pxy     <- qs::qread(d$file_coord)$states
        # Unwrap coast
        # * This is very fast (2 milliseconds)
        # * So there is no need to implement chunking
        coast_s <- terra::unwrap(coast_s_w)
        # Compute area of MCP by timestep & return data.table (~1.5 mins)
        cl_lapply(split(pxy, pxy$timestep), function(pxy_for_t) {
          # pxy_for_t <- pxy[timestep == 3L, ]
          if (FALSE) {
            terra::plot(map)
            points(pxy_for_t$x, pxy_for_t$y)
          }
          coord_n <- collapse::fnunique(rleid(pxy_for_t$x, pxy_for_t$y))
          if (coord_n == 1L) {
            coord_km2 <- 0
          } else {
            coord_km2 <- 
              cbind(pxy_for_t$x, pxy_for_t$y) |> 
              terra::vect(crs = terra::crs(coast_s)) |> 
              terra::convHull() |> 
              terra::erase(coast_s) |>
              terra::expanse(unit = "km")
          }
          data.table(file_coord = d$file_coord, 
                     timestep   = pxy_for_t$timestep[1], 
                     coord_n    = coord_n,
                     coord_km2  = coord_km2
          )
        }) |> rbindlist()
      }
    }) |> rbindlist()
  toc()
  qs::qsave(particle_mcps, here_data("output", "analysis-summary", "particle-mcps.qs"))
  
}


###########################
#### Compute residency (~22 s)

if (TRUE) {
  
  iteration_res <- copy(iteration)
  iteration_res[, file := file.path(iteration$folder_ud, "spatstat", "h", "ud.tif")]
  iteration_res[, file_exists := file.exists(file)]
  residency <- lapply_estimate_residency_ud(files = iteration_res$file[iteration_res$file_exists])
  # Write output
  residency <- 
    left_join(iteration_res, residency, by = "file") |> 
    mutate(algorithm = dataset) |>
    select(individual_id, month_id, unit_id, algorithm, sensitivity, zone, time) |> 
    arrange(individual_id, month_id, unit_id, algorithm, sensitivity, zone) |>
    as.data.table()
  qs::qsave(residency, here_data("output", "analysis-summary", "residency-patter.qs"))
  
}


#### End of code.
###########################
###########################