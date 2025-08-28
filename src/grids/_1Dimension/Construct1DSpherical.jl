




struct SphericalGrid1D
    coord1::UniformAxis
    bounds::Tuple{Float64, Float64}
    ghost_bounds::Tuple{Float64, Float64}
    units::String
end

function Construct1DSpherical(domain_length::Float64, zones::Int, ghost_zones::Int, grid_center::Float64, units::String = "cm")
    coord1 = ConstructUniformAxis(domain_length, zones, ghost_zones, grid_center, :spherical)
    bounds = (0.0 + coord1.spacing/2 + grid_center, domain_length + grid_center)
    ghost_bounds = (coord1.spacing/2 + grid_center - ghost_zones * coord1.spacing,
                    domain_length + grid_center + ghost_zones * coord1.spacing)
    return SphericalGrid1D(coord1, bounds, ghost_bounds, units)
end