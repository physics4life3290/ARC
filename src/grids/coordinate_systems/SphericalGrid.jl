



struct SphericalGrid{N,T<:AbstractFloat}
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

function ConstructSphericalGrid(
        domain_lengths::NTuple{N,T},   # (r_max,) or (r_max, θ_max) or (r_max, θ_max, φ_max)
        zones::NTuple{N,Int},
        ghost_zones::Int,
        grid_center::NTuple{N,T}
    ) where {N,T<:AbstractFloat}

    # Choose axis types dynamically
    axis_types = if N == 1
        (:radial,)
    elseif N == 2
        (:radial, :angular)
    else
        (:radial, :angular, :angular)
    end

    # 1. Build coordinate axes
    Axes = ntuple(i -> ConstructUniformAxis(
        domain_lengths[i],
        zones[i],
        ghost_zones,
        grid_center[i],
        axis_types[i]
    ), N)

    # 2. Compute physical bounds
    Bounds = ntuple(i -> begin
        at = axis_types[i]
        if at == :radial
            (grid_center[i], grid_center[i] + domain_lengths[i])  # r ≥ 0
        elseif at == :polar
            (zero(T), domain_lengths[i])                          # θ in [0, θ_max]
        else
            (zero(T), domain_lengths[i])                          # φ in [0, φ_max]
        end
    end, N)

    # 3. Domain = same as Bounds
    Domain = Bounds

    # 4. Ghost bounds using spacing
    Ghost_Bounds = ntuple(i -> (
        Bounds[i][1] - Axes[i].spacing * ghost_zones,
        Bounds[i][2] + Axes[i].spacing * ghost_zones
    ), N)

    # 5. Total zones including ghost cells
    Total_Zones = ntuple(i -> zones[i] + 2*ghost_zones, N)

    return SphericalGrid{N,T}(
        Axes,           # axes
        Domain,         # domain
        Bounds,         # physical bounds
        Ghost_Bounds,   # ghost bounds
        zones,          # zones
        ghost_zones,    # ghost zones
        Total_Zones,    # total zones
        grid_center,    # origin
        N               # dimension
    )
end
