struct ModelObsCaptureContainer <: Patter.ModelObs
    sensor_id::Int64
    capture_x::Float64
    capture_y::Float64
    radius::Float64
end

function Patter.logpdf_obs(state::State, model::ModelObsCaptureContainer, t::Int64, obs::Int64)
    # Calculate distance between particle (state) and capture location
    dist = Patter.distance(state.x, state.y, model.capture_x, model.capture_y)
    # Only particles within model.radius are permitted
    return ifelse(dist <= model.radius, 0.0, -Inf)
end
