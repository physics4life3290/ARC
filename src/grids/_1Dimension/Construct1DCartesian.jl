




struct CartesianGrid1D
    xcoord::UniformAxis
    bounds::Tuple{Float64, Float64}
    ghost_bounds::Tuple{Float64, Float64}
    units::String
end


function Construct1DCartesian(domain_length::Float64, zones::Int, ghost_zones::Int, grid_center::Float64, units::String = "cm")
    xcoord = ConstructUniformAxis(domain_length, zones, ghost_zones, grid_center, :cartesian)
    bounds = (-domain_length/2 + grid_center, domain_length/2 + grid_center)
    ghost_bounds = (-domain_length/2 + grid_center - ghost_zones * xcoord.spacing, 
                    domain_length/2 + grid_center + ghost_zones * xcoord.spacing)
    return CartesianGrid1D(xcoord, bounds, ghost_bounds, units)
end