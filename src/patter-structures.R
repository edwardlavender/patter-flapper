#' ModelMoveFlapperRW

patter_ModelMoveRW <- function(sim) {
  stopifnot(all(c("k1", "theta1", "k2", "theta2", "mobility") %in% colnames(sim)))
  move_flapper_rw(mobility = sim$mobility, 
                  dbn_length_rest = glue::glue("truncated(Cauchy({sim$k1}, {sim$theta1}), lower = 0.0, upper = {sim$mobility})"),
                  dbn_length_active = glue::glue("truncated(Cauchy({sim$k2}, {sim$theta2}), lower = 0.0, upper = {sim$mobility})"),
                  dbn_angle = "Uniform(-pi, pi)")
}

move_flapper_rw <- function(mobility = "1095.0", 
                            dbn_length_rest = "truncated(Cauchy(0.0, 5.0), lower = 0.0, upper = 1095.0)", 
                            dbn_length_active = "truncated(Cauchy(5.0, 100.0), lower = 0.0, upper = 1095.0)", 
                            dbn_angle = "Uniform(-pi, pi)") {
  julia_check_exists("env")
  glue("ModelMoveFlapperRW(env, {mobility}, {dbn_length_rest}, {dbn_length_active}, {dbn_angle});")
}

# Check behaviour of Cauchy distribution matches behaviour in R: yes. 
# using Plots
# using Distributions 
# dist = truncated(Cauchy(0.0, 5.0), lower = 0.0, upper = 1095.0)
# histogram(rand(dist, 100000), bins = 100)
# dist = truncated(Cauchy(5.0, 100.0), lower = 0.0, upper = 1095.0)
# histogram(rand(dist, 100000), bins = 100)

ModelMoveFlapperRW <- 
  '
  struct ModelMoveFlapperRW{T, U, V, W, X} <: Patter.ModelMove
    map::T
    mobility::U
    dbn_length_rest::V
    dbn_length_active::W
    dbn_angle::X
  end
  '

#' simulate_step() method
# This function assumes that the 'behaviour' vector (1, 2) exists in Julia when sourced
simulate_step.ModelMoveFlapperRW  <- 
  '
  function Patter.simulate_step(state::StateXY, model_move::ModelMoveFlapperRW, t::Int64, behaviour::Vector{Int64} = behaviour)
  
    # Simulate step length, depending on behaviour
    if behaviour[t] == 1
      length = rand(model_move.dbn_length_rest)
    elseif behaviour[t] == 2
      length = rand(model_move.dbn_length_active)
    else
      error("behaviour vector must be 1/2.")
    end 
    
    angle     = rand(model_move.dbn_angle)
    x         = state.x + length * cos(angle)
    y         = state.y + length * sin(angle)
    map_value = Patter.extract(model_move.map, x, y)
    StateXY(map_value, x, y)
  end 
'

#' logpdf_step() method
# This function assumes that the 'behaviour' vector (1, 2) exists in Julia when sourced

logpdf_step.ModelMoveFlapperRW  <- 
  '
  function Patter.logpdf_step(state_from::StateXY, state_to::StateXY, model_move::ModelMoveFlapperRW, 
                       t::Int64, length::Float64, angle::Float64, behaviour::Vector{Int64} = behaviour) 
    if behaviour[t] == 1
      length_logpdf = logpdf(model_move.dbn_length_rest, length)
    elseif behaviour[t] == 2
      length_logpdf = logpdf(model_move.dbn_length_active, length)
    end 
    length_logpdf + logpdf(model_move.dbn_angle, angle)
  end 
  '


#' ModelMoveFlapperCRW 
#' (Used to test whether a CRW improves convergence properties)

patter_ModelMoveCRW <- function(sim) {
  stopifnot(all(c("k1", "theta1", "k2", "theta2", "mobility") %in% colnames(sim)))
  move_flapper_crw(mobility = sim$mobility, 
                  dbn_length_rest = glue::glue("truncated(Cauchy({sim$k1}, {sim$theta1}), lower = 0.0, upper = {sim$mobility})"),
                  dbn_length_active = glue::glue("truncated(Cauchy({sim$k2}, {sim$theta2}), lower = 0.0, upper = {sim$mobility})"))
}

rho <- "1.5"
move_flapper_crw <- function(mobility = "1095.0", 
                             dbn_length_rest = "truncated(Cauchy(0.0, 5.0), lower = 0.0, upper = 1095.0)", 
                             dbn_length_active = "truncated(Cauchy(5.0, 100.0), lower = 0.0, upper = 1095.0)", 
                             dbn_angle_delta = glue("Normal(0.0, {rho})")) {
  julia_check_exists("env")
  glue("ModelMoveFlapperCRW(env, {mobility}, {dbn_length_rest}, {dbn_length_active}, {dbn_angle_delta});")
}

StateXYD <- 
  '
struct StateXYD <: Patter.State
    map_value::Float64
    x::Float64
    y::Float64
    angle::Float64  
end
  '

states_init.StateXYD <- 
  '
  function Patter.states_init(state_type::Type{StateXYD}, coords)
    N = nrow(coords)
    coords.angle = rand(N) .* 2 .* π
    return coords
  end 
  '

ModelMoveFlapperCRW <- 
  '
  struct ModelMoveFlapperCRW{T, U, V, W, X} <: Patter.ModelMove
    map::T
    mobility::U
    dbn_length_rest::V
    dbn_length_active::W
    dbn_angle_delta::X
  end
  '

simulate_step.ModelMoveFlapperCRW <- 
  '
  function Patter.simulate_step(state::StateXYD, model_move::ModelMoveFlapperCRW, t::Int64, behaviour::Vector{Int64} = behaviour)
  
    # Simulate step length, depending on behaviour
    if behaviour[t] == 1
      length = rand(model_move.dbn_length_rest)
    elseif behaviour[t] == 2
      length = rand(model_move.dbn_length_active)
    else
      error("behaviour vector must be 1/2.")
    end 
    
    angle     = state.angle + rand(model_move.dbn_angle_delta)
    x         = state.x + length * cos(angle)
    y         = state.y + length * sin(angle)
    map_value = Patter.extract(model_move.map, x, y)
    StateXYD(map_value, x, y, angle)
  end 
'

# TO DO
logpdf_step.ModelMoveFlapperCRW <- 
  '
  function Patter.logpdf_step(state_from::StateXYD, state_to::StateXYD, model_move::ModelMoveFlapperCRW, 
                              t::Int64, length::Float64, angle::Float64, behaviour::Vector{Int64} = behaviour) 
                              
    # Compute logpdf_length                          
    if behaviour[t] == 1
      length_logpdf = logpdf(model_move.dbn_length_rest, length)
    elseif behaviour[t] == 2
      length_logpdf = logpdf(model_move.dbn_length_active, length)
    end 
    
    # Compute logpdf_angle
    angle_delta = Patter.abs_angle_difference(angle, state_from.angle)
    angle_logpdf = logpdf(model_move.dbn_angle_delta, angle_delta)
    
    # Return (unnnormalised) density of length & turning angle
    return length_logpdf + angle_logpdf
    
  end 
  '

#' ModelMoveFlapper (RW or CRW)
state_flapper                    <- "StateXYD"
patter_ModelMove                 <- patter_ModelMoveCRW
move_flapper                     <- move_flapper_crw
ModelMoveFlapper                 <- ModelMoveFlapperCRW
simulate_step.ModelMoveFlapper   <- simulate_step.ModelMoveFlapperCRW
logpdf_step.ModelMoveFlapper     <- logpdf_step.ModelMoveFlapperCRW

#' Helper routine to set initial model components 

set_model_move_components <- function() {
  # (optional) State State methods for StateXYD (for CRW)
  JuliaCall::julia_command(StateXYD)
  JuliaCall::julia_command(states_init.StateXYD)
  # Set movement model (RW or CRW)
  JuliaCall::julia_command(ModelMoveFlapper)
  nothing()
}

update_model_move_components <- function() {
  JuliaCall::julia_command(simulate_step.ModelMoveFlapper)
  JuliaCall::julia_command(logpdf_step.ModelMoveFlapper)
  nothing()
}

#' Helper routine to determine if RW or CRW

model_move_is_crw <- function() {
  "dbn_angle_delta" %in% names(formals(move_flapper))
}
