




struct CylindricalGrid{N,T<:AbstractFloat}
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


function ConstructCylindricalGrid(
        domain_lengths::NTuple{N,T},   # e.g. (r_max,) or (r_max, phi_max) or (r_max, phi_max, z_len)
        zones::NTuple{N,Int},
        ghost_zones::Int,
        grid_center::NTuple{N,T}
    ) where {N,T<:AbstractFloat}

    # choose axis types depending on dimensionality
    axis_types = if N == 1
        (:radial,)
    elseif N == 2
        (:radial, :angular)
    else
        (:radial, :angular, :linear)
    end

    # build axes (ConstructUniformAxis expects (domain, zones, ghost_zones, center, axis_type))
    Axes = ntuple(i -> ConstructUniformAxis(domain_lengths[i], zones[i], ghost_zones, grid_center[i], axis_types[i]), N)

    # compute physical bounds consistent with axis type and ConstructUniformAxis behaviour
    Bounds = ntuple(i -> begin
            at = axis_types[i]
            if at == :radial
                (grid_center[i], grid_center[i] + domain_lengths[i])        # r: [center, center + rmax]
            elseif at == :angular
                (zero(T), domain_lengths[i])                                # φ: [0, φ_max]
            else # :linear
                (grid_center[i] - domain_lengths[i]/2, grid_center[i] + domain_lengths[i]/2)
            end
        end, N)

    # domain field in struct expects (min,max) tuples for each axis — use same as Bounds
    Domain = Bounds

    # ghost bounds using axis spacing from each UniformAxis
    Ghost_Bounds = ntuple(i -> (Bounds[i][1] - Axes[i].spacing * ghost_zones,
                                Bounds[i][2] + Axes[i].spacing * ghost_zones), N)

    # total zones including ghost cells
    Total_Zones = ntuple(i -> zones[i] + 2*ghost_zones, N)

    return CylindricalGrid{N,T}(
        Axes,          # axes
        Domain,        # domain (min,max)
        Bounds,        # bounds
        Ghost_Bounds,  # ghost_bounds
        zones,         # zones
        ghost_zones,   # ghost_zones
        Total_Zones,   # total_zones
        grid_center,   # origin
        N              # dimension
    )
end

