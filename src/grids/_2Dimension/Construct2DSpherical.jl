









struct SphericalGrid2D
    rcoord::UniformAxis
    thetacoord::UniformAxis
    rbounds::Tuple{Float64, Float64}
    thetabounds::Tuple{Float64, Float64}
    rghost_bounds::Tuple{Float64, Float64}
    thetaghost_bounds::Tuple{Float64, Float64}
    units::String
end

function Construct2DSpherical(domain_length::Float64, zones::Int, θzones::Int, ghost_zones::Int, grid_center::Float64, units::String = "cm")
    rcoord = ConstructUniformAxis(domain_length, zones, ghost_zones, grid_center, :spherical)
    thetacoord = ConstructUniformAxis(Float64(π), θzones, ghost_zones, grid_center, :cartesian)
    rbounds = (0.0 + rcoord.spacing/2, domain_length - rcoord.spacing/2)
    thetabounds = (0.0 + thetacoord.spacing/2, π - thetacoord.spacing/2)
    rghost_bounds = (-domain_length/2 + grid_center - ghost_zones * rcoord.spacing,
                    domain_length/2 + grid_center + ghost_zones * rcoord.spacing)
    thetaghost_bounds = (-π/2 - ghost_zones * thetacoord.spacing,
                    π/2 + ghost_zones * thetacoord.spacing)
    return SphericalGrid2D(rcoord, thetacoord, rbounds, thetabounds, rghost_bounds, thetaghost_bounds, units)
end