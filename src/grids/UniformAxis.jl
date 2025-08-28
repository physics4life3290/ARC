




struct UniformAxis
    centers::Vector{Float64}
    interfaces::Vector{Float64}
    all_centers::Vector{Float64}
    all_interfaces::Vector{Float64}
    spacing::Float64
    total_zones::Int
    zones::Int
    ghost_zones::Int
    coord_system::Symbol
end

function ConstructUniformAxis(domain_length::Float64, zones::Int, ghost_zones::Int, grid_center::Float64, coord_system::Symbol)
    total_zones = zones + 2 * ghost_zones
    spacing = domain_length / zones

    if coord_system == :cartesian
        centers = LinRange(-domain_length/2 + grid_center, domain_length/2 + grid_center, zones)
        interfaces = LinRange(-domain_length/2 + grid_center - spacing/2, domain_length/2 + grid_center + spacing/2, zones + 1)

        all_centers = LinRange(-domain_length/2 + grid_center - ghost_zones * spacing,
                               domain_length/2 + grid_center + ghost_zones * spacing,
                               total_zones)

        all_interfaces = LinRange(-domain_length/2 + grid_center - (ghost_zones + 0.5) * spacing,
                                  domain_length/2 + grid_center + (ghost_zones + 0.5) * spacing,
                                  total_zones + 1)

    elseif coord_system == :spherical || coord_system == :cylindrical

        r_min = 1E-12
        r_max = domain_length

        centers = LinRange(r_min + grid_center, r_max + grid_center, zones)
        interfaces = LinRange(r_min + grid_center - spacing/2, r_max + grid_center + spacing/2, zones + 1)

        all_centers = LinRange(r_min + grid_center - ghost_zones * spacing,
                               r_max + grid_center + ghost_zones * spacing,
                               total_zones)

        all_interfaces = LinRange(r_min + grid_center - (ghost_zones + 0.5) * spacing,
                                  r_max + grid_center + (ghost_zones + 0.5) * spacing,
                                  total_zones + 1)
    else
        error("Unsupported coordinate system: $coord_system")
    end

    return UniformAxis(centers, interfaces, all_centers, all_interfaces, spacing, total_zones, zones, ghost_zones, coord_system)
end
