




# ----------------------
# Abstract Grid
# ----------------------
abstract type _Grid end

# ----------------------
# Axes
# ----------------------
struct UniformAxis{T}
    centers::Vector{T}
    faces::Vector{T}
    dx::T
    zones::Int32
    bounds::Tuple{T,T}
end

# outer constructor: automatically converts ranges to vectors
UniformAxis(origin::T, dx::T, zones::Int32) where {T} = UniformAxis(
    collect(origin .+ dx*(0.5:(zones-0.5))),
    collect(origin .+ dx*(0:zones)),
    dx,
    zones,
    (origin, origin + dx*zones)
)


struct NonUniformAxis{T}
    centers::Vector{T}
    faces::Vector{T}
    dx::Vector{T}
    zones::Int32
    bounds::Tuple{T,T}
end

struct AdaptiveAxis{T}
    centers::Vector{T}
    faces::Vector{T}
    dx::Vector{T}
    zones::Int32
    bounds::Tuple{T,T}
end

# ----------------------
# Grids
# ----------------------
struct UniformGrid{T,N} <: _Grid
    axes::NTuple{N, UniformAxis{T}}
    origin::NTuple{N,T}
    ghost_zones::Int32
    total_zones::NTuple{N, Int32}
    active_zones::NTuple{N, UnitRange{Int32}}
    domain::NTuple{N, T}
end

struct NonUniformGrid{T,N} <: _Grid
    axes::NTuple{N, NonUniformAxis{T}}
    origin::NTuple{N,T}
    ghost_zones::Int32
    total_zones::NTuple{N,Int32}
    active_zones::NTuple{N,UnitRange{Int32}}
    domain::NTuple{N, T}
end

mutable struct Block{T}
    data::Vector{Vector{T}}
    dx::T
    level::Int32
    children::Vector{Block{T}}
    x0::T
    nx_internal::Int32
end

struct AdaptiveGrid{T,N} <: _Grid
    root_blocks::Vector{Block{T}}
    origin::NTuple{N,T}
    domain::NTuple{N,T}
    ghost_zones::Int32
    max_level::Int32
end