




struct SphericalGrid1D
    rcoord::UniformAxis
    bounds::Tuple{Float64, Float64}
    ghost_bounds::Tuple{Float64, Float64}
    units::String
end

function Construct1DSpherical(domain_length::Float64, zones::Int, ghost_zones::Int, grid_center::Float64, units::String = "cm")
    rcoord = ConstructUniformAxis(domain_length, zones, ghost_zones, grid_center, :spherical)
    bounds = (0.0 + rcoord.spacing/2, domain_length - rcoord.spacing/2)
    ghost_bounds = (-domain_length/2 + grid_center - ghost_zones * rcoord.spacing,
                    domain_length/2 + grid_center + ghost_zones * rcoord.spacing)
    return SphericalGrid1D(rcoord, bounds, ghost_bounds, units)
end