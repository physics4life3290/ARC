




struct CartesianGrid1D{T<:AbstractFloat}
    coord1::UniformAxis   # stays concrete
    bounds::Tuple{T,T}
    ghost_bounds::Tuple{T,T}
    units::String
end

@inline function Construct1DCartesian(domain_length::T, zones::Int, ghost_zones::Int,
                                      grid_center::T, units::String = "cm") where {T<:AbstractFloat}
    coord1 = ConstructUniformAxis(domain_length, zones, ghost_zones, grid_center, :cartesian)

    half    = domain_length / 2
    lo      = grid_center - half
    hi      = grid_center + half
    spacing = coord1.spacing

    bounds       = (lo, hi)
    ghost_bounds = (lo - ghost_zones * spacing,
                    hi + ghost_zones * spacing)

    return CartesianGrid1D{T}(coord1, bounds, ghost_bounds, units)
end
