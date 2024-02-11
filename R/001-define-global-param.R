###########################
###########################
#### define-global-param.R

#### Aims
# 1) Define parameters

#### Prerequisites
# 1) flapper_appl project
# 2) These parameters are informed by prepare-patter-param.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(plotly)
library(patter)


###########################
###########################
#### Define parameters

#### Step duration
step <- "2 mins"

#### Detection range
detection_range <- 1000

#### Mobility 
mobility <- 750

#### Truncated gamma movement model
sh <- 1
sc <- 250
hist(rtruncgamma(1e5, .shape = sh, .scale = sc, .mobility = mobility), 
     probability = TRUE, breaks = 100, xlim = c(0, mobility + 200))
x <- seq(0, mobility, by = 1)
dtruncgamma(0, .shape = sh, .scale = sc, .mobility = mobility)
dtruncgamma(1, .shape = sh, .scale = sc, .mobility = mobility)
y <- dtruncgamma(x, .shape = sh, .scale = sc, .mobility = mobility)
data <- data.frame(x = x, 
                   y = y)
plot_ly(data = data, x = ~x, y = ~y) |> 
  add_lines() |>
  layout(yaxis = list(range = c(0, 0.005)))

#### Uniform movement model 
hist(runif(1e5, 0, mobility))

### (original) Depth-error function 
calc_depth_error <- function(depth) {
  e <- 4.77 + 2.5 + sqrt(0.5^2 + (0.013 * depth)^2)
  matrix(c(-(e + 5), e), nrow = 2)
}
calc_depth_error <- Vectorize(calc_depth_error)

#### Collate pars
pars <- 
  list(
    flapper = list(
      step = step,
      detection_range = 750, 
      mobility = 500, 
      calc_depth_error = calc_depth_error), 
    patter = list(
      step = step,
      detection_range = detection_range, 
      shape = sh, 
      scale = sc, 
      mobility = mobility, 
      # Define the model we will use
      # * The "truncated gamma" or "uniform" model above
      # * We include all parameters in this list either way
      # * (to facilitate exploration of alternative choices)
      model = "uniform"
      )
  )

#### Save pars
saveRDS(pars, here_data("input", "pars.rds"))


#### End of code. 
###########################
###########################