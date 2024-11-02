#' ModelMoveFlapper

move_flapper <- function(mobility = "1095.0", 
                         dbn_length_rest = "truncated(Cauchy(0.0, 5.0), lower = 0.0, upper = 1095.0)", 
                         dbn_length_active = "truncated(Cauchy(5.0, 100.0), lower = 0.0, upper = 1095.0)", 
                         dbn_angle = "Uniform(-pi, pi)") {
  julia_check_exists("env")
  glue("ModelMoveFlapper(env, {mobility}, {dbn_length_rest}, {dbn_length_active}, {dbn_angle});")
}

# Check behaviour of Cauchy distribution matches behaviour in R: yes. 
# using Plots
# using Distributions 
# dist = truncated(Cauchy(0.0, 5.0), lower = 0.0, upper = 1095.0)
# histogram(rand(dist, 100000), bins = 100)
# dist = truncated(Cauchy(5.0, 100.0), lower = 0.0, upper = 1095.0)
# histogram(rand(dist, 100000), bins = 100)

ModelMoveFlapper <- 
  '
  struct ModelMoveFlapper{T, U, V, W, X} <: Patter.ModelMove
    map::T
    mobility::U
    dbn_length_rest::V
    dbn_length_active::W
    dbn_angle::X
  end
  '

#' simulate_step() method
# This function assumes that the 'behaviour' vector (1, 2) exists in Julia when sourced
simulate_step.ModelMoveFlapper <- 
  '
  function Patter.simulate_step(state::StateXY, model_move::ModelMoveFlapper, t::Int64, behaviour::Vector{Int64} = behaviour)
  
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

logpdf_step.ModelMoveFlapper <- 
  '
  function Patter.logpdf_step(state_from::StateXY, state_to::StateXY, model_move::ModelMoveFlapper, 
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

move_flapper_crw <- function(mobility = "1095.0", 
                             dbn_length_rest = "truncated(Cauchy(0.0, 5.0), lower = 0.0, upper = 1095.0)", 
                             dbn_length_active = "truncated(Cauchy(5.0, 100.0), lower = 0.0, upper = 1095.0)", 
                             dbn_angle_delta = "Normal(0.0, 1.2)") {
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

#' simulate_step() method
# This function assumes that the 'behaviour' vector (1, 2) exists in Julia when sourced
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