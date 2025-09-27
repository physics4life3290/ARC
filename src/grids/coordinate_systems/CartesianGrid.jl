




struct CartesianGrid{N,T<:AbstractFloat}
    axes::NTuple{N, UniformAxis{T}}
    domain::NTuple{N, Tuple{T,T}}
    bounds::NTuple{N, Tuple{T,T}}
    ghost_bounds::NTuple{N, Tuple{T,T}}
    zones::NTuple{N, Int}
    ghost_zones::Int
    total_zones::NTuple{N, Int}
    origin::NTuple{N, T}
    dimension::Int
end

function ConstructCartesianGrid(
        domain_lengths::NTuple{N,T},
        zones::NTuple{N,Int},
        ghost_zones::Int,
        grid_center::NTuple{N,T}
    ) :: CartesianGrid{N,T} where {N,T<:AbstractFloat}

    # 1. Build coordinate axes
    Axes = ntuple(i -> ConstructUniformAxis(domain_lengths[i], zones[i], ghost_zones, grid_center[i], :linear), N)

    # 2. Physical bounds
    Bounds = ntuple(i -> (grid_center[i] - domain_lengths[i]/2,
                          grid_center[i] + domain_lengths[i]/2), N)

    # 3. Ghost bounds using axis spacing
    Ghost_Bounds = ntuple(i -> (Bounds[i][1] - Axes[i].spacing * ghost_zones,
                                Bounds[i][2] + Axes[i].spacing * ghost_zones), N)

    # 4. Total zones including ghost cells
    Total_Zones = ntuple(i -> zones[i] + 2*ghost_zones, N)

    return CartesianGrid{N,T}(Axes, Bounds, Bounds, Ghost_Bounds, zones, ghost_zones, Total_Zones, grid_center, N)
end
