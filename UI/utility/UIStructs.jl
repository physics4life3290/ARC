




struct PrimaryInput
    problem::Symbol
    mode::Symbol
    dimension::Int
    coordinate_system::Symbol
    solver::Symbol
    boundary_condition::Symbol
    filename::String
end

# We will build this struct as we go up in dimensions
# We can use a constructor that takes in primaryinput
# this will build a multiD grid each with fully 
# customizable fields.
struct GridInput
    domain::Float64
    grid_center::Float64
    zones::Int
    ghost_zones::Int
    total_zones::Int
    coord_min::Float64
    coord_max::Float64
end

struct SolverInput
    cfl::Float64
    t_final::Float64
    flattening::Bool
    steepening::Bool
    limiter::Union{Symbol, Nothing}
    reconstruction::Union{Symbol, Nothing}
    riemanntype::Union{Symbol, Nothing}
end

struct ShockTubeInput
    wall_positions::Vector{Float64}
    states::Vector{PrimitiveVariables}
    gamma::Float64
end

struct BlastRegion
    center::Float64
    radius::Float64
    energy::Float64
end

struct BlastParameters
    states::Vector{PrimitiveVariables}
    blasts::Vector{BlastRegion}
end

struct BlastWaveInput
    parameters::BlastParameters
    wall_positions::Vector{Float64}
    ambient_states::Vector{PrimitiveVariables}
    gamma::Float64
end



struct UserInput
    Primary_Input::PrimaryInput
    Grid_Input::GridInput
    Solver_Input::SolverInput
    Secondary_Input::Union{ShockTubeInput, BlastWaveInput}
end