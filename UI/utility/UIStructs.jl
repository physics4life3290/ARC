




struct PrimaryInput
    problem::Symbol
    mode::Symbol
    features::Vector{Symbol}
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
struct GridInput{N}
    domain::NTuple{N,Float64}
    grid_center::NTuple{N,Float64}
    zones::NTuple{N,Int}
    ghost_zones::Int
    total_zones::NTuple{N,Int}
    coord_min::NTuple{N,Float64}
    coord_max::NTuple{N,Float64}
end


struct SolverInput
    cfl::Float64
    t_final::Float64
    split_choice::Symbol
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