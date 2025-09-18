




# Shared base state for all levels
struct BaseState{T<:AbstractVector}
    ρ::T
    p::T
    γ::Float64
    grid_points::T
    coordinate_system::Symbol
end

# ---------------------
# Progressive solver access levels
# ---------------------

struct LevelI{T<:AbstractVector}
    state::BaseState{T}
end

struct LevelII{T<:AbstractVector}
    state::BaseState{T}
    solver::Symbol
    Cfl::Float64
    boundary_condition::Symbol
end

struct LevelIII{T<:AbstractVector}
    state::BaseState{T}
    solver::Symbol
    Cfl::Float64
    boundary_condition::Symbol
    mode::Symbol
    features::Vector{Symbol}
end

struct LevelIV{T<:AbstractVector}
    state::BaseState{T}
    solver::Symbol
    Cfl::Float64
    boundary_condition::Symbol
    mode::Symbol
    features::Vector{Symbol}
    reconstruction::Symbol
    steepening::Bool
    flattening::Bool
    limiter::Symbol
    riemann_solver::Symbol
end

