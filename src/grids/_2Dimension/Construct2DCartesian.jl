




struct CartesianGrid2D
    xcoord::UniformAxis
    ycoord::UniformAxis
    xbounds::Tuple{Float64, Float64}
    ybounds::Tuple{Float64, Float64}
    xghost_bounds::Tuple{Float64, Float64}
    yghost_bounds::Tuple{Float64, Float64}
    grid_center::Float64
    cell_areas::Matrix{Float64}
    units::String
end


function Construct2DCartesian(
    xdomain_length::Float64,
    ydomain_length::Float64,
    xzones::Int,
    yzones::Int,
    ghost_zones::Int,
    grid_center::Float64,
    units::String
)
    xcoord = ConstructUniformAxis(xdomain_length, xzones, ghost_zones, grid_center, :cartesian)
    ycoord = ConstructUniformAxis(ydomain_length, yzones, ghost_zones, grid_center, :cartesian)

    xbounds = (-xdomain_length/2 + grid_center, xdomain_length/2 + grid_center)
    ybounds = (-ydomain_length/2 + grid_center, ydomain_length/2 + grid_center)

    cell_areas = xcoord.spacing * ycoord.spacing .* ones(xzones, yzones)

    xghost_bounds = (-xdomain_length/2 + grid_center - ghost_zones * xcoord.spacing, 
                     xdomain_length/2 + grid_center + ghost_zones * xcoord.spacing)

    yghost_bounds = (-ydomain_length/2 + grid_center - ghost_zones * ycoord.spacing, 
                     ydomain_length/2 + grid_center + ghost_zones * ycoord.spacing)

    return CartesianGrid2D(xcoord, ycoord, xbounds, ybounds, xghost_bounds, yghost_bounds, grid_center, cell_areas, units) 
end

