




struct UniformAxis{T}
    PhysicalCenters::Vector{T}
    PhysicalInterfaces::Vector{T}
    SimulationCenters::Vector{T}
    SimulationInterfaces::Vector{T}
    spacing::T
    AxisCenter::T
end

function ConstructUniformAxis(domain::T, zones::Int, ghost_zones::Int, grid_center::T, axis_type::Symbol) where T
    
    spacing = domain / zones
    total_zones = zones + 2 * ghost_zones

    if axis_type == :linear
        
        half_domain = domain / 2
        min_physical = -half_domain + grid_center
        max_physical = half_domain + grid_center

        # Physical centers and interfaces
        physical_centers = collect(range(min_physical, max_physical, length=zones))
        physical_interfaces = collect(range(min_physical + spacing/2, max_physical - spacing/2, length=zones-1))

        # Ghost-extended simulation bounds
        min_ghost_bound = first(physical_centers) - spacing * ghost_zones
        max_ghost_bound = last(physical_centers) + spacing * ghost_zones

        simulation_centers = collect(range(min_ghost_bound, max_ghost_bound, length=total_zones))
        simulation_interfaces = collect(range(min_ghost_bound + spacing/2, max_ghost_bound - spacing/2, length=total_zones-1))

        return UniformAxis{T}(physical_centers, physical_interfaces, simulation_centers, simulation_interfaces, spacing, grid_center)
    
    elseif axis_type == :radial

        min_physical = grid_center
        max_physical = min_physical + domain
        half_domain = (max_physical - min_physical) / 2

        physical_centers = collect(range(min_physical, max_physical, length=zones))
        physical_interfaces = collect(range(min_physical+spacing/2, max_physical-spacing/2, length=zones-1))

        total_zones = zones + 2 * ghost_zones
        min_ghost_bound = first(physical_centers) - spacing * ghost_zones
        max_ghost_bound = last(physical_centers) + spacing * ghost_zones

        simulation_centers = collect(range(min_ghost_bound, max_ghost_bound, length=total_zones))
        simulation_interfaces = collect(range(min_ghost_bound+spacing/2, max_ghost_bound-spacing/2, length=total_zones-1))

        return UniformAxis{T}(physical_centers, physical_interfaces, simulation_centers, simulation_interfaces, spacing, grid_center)

    elseif axis_type == :angular

        min_physical = 0.0
        max_physical = 3.141592653589793
        half_domain = (max_physical - min_physical) / 2

        physical_centers = collect(range(min_physical, max_physical, length=zones))
        physical_interfaces = collect(range(min_physical+spacing/2, max_physical-spacing/2, length=zones-1))

        min_ghost_bound = first(physical_centers) - spacing * ghost_zones
        max_ghost_bound = last(physical_centers) + spacing * ghost_zones

        total_zones = zones + 2 * ghost_zones

        simulation_centers = collect(range(min_ghost_bound, max_ghost_bound, length=total_zones))
        simulation_interfaces = collect(range(min_ghost_bound+spacing/2, max_ghost_bound-spacing/2, length=total_zones-1))
        return UniformAxis{T}(physical_centers, physical_interfaces, simulation_centers, simulation_interfaces, spacing, grid_center)
    end

end


#=
zones = 10
ghost_zones = 3
domain = 2.0
name = "x-Axis"
grid_center = 0.0

# Grid takes in domain::NTuple(N, T), zones::NTuple{N,T}, ghost_zones, grid_center::NTuple{N,T}, names::NTuple{N,String}

coord = ConstructUniformAxis(domain, zones, ghost_zones, grid_center, name)

println(coord.PhysicalCenters)
println(coord.PhysicalInterfaces)
println(coord.SimulationCenters)
println(coord.SimulationInterfaces)
println(coord.spacing)
println(coord.AxisCenter)
println(coord.name)

struct Grid{N, T}
    axes::NTuple{N, UniformAxis}
    domain::NTuple{N, T}
    bounds::NTuple{N, Tuple{T, T}}
    ghost_bounds::NTuple{N, Tuple{T, T}}
    zones::NTuple{N, Int}
    ghost_zones::Int
    total_zones::Int
    origin::NTuple{N, T}
    dimension::N
    coordinate_system::Symbol
end
=#