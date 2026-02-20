"""
Uniform 1D grid axis with ghost zones.

Fields:
- centers: cell centers including ghost zones
- faces: cell interfaces including ghost zones
- dx: uniform spacing
- interior_zones: number of physical (non-ghost) cells
- ghost_zones: ghost zones per side
- ilo: first interior index
- ihi: last interior index
- bounds: (xmin, xmax) of physical domain
"""
struct UniformAxis{T}
    centers::Vector{T}
    faces::Vector{T}
    dx::T
    interior_zones::Int
    ghost_zones::Int
    ilo::Int
    ihi::Int
    bounds::Tuple{T,T}
end


function UniformAxis(bounds::Tuple{T,T}, zones::Int; ghost_zones::Int=3) where T

    @assert zones > 0 "zones must be positive"
    @assert ghost_zones >= 0 "ghost_zones must be non-negative"

    xmin, xmax = bounds
    @assert xmax > xmin "invalid bounds"

    dx = (xmax - xmin) / zones

    total_zones = zones + 2*ghost_zones

    centers = Vector{T}(undef, total_zones)
    faces   = Vector{T}(undef, total_zones + 1)

    # First interior center
    first_center = xmin + dx/2

    # First center including ghosts
    start_center = first_center - ghost_zones*dx

    # Fill centers
    @inbounds for i in 1:total_zones
        centers[i] = start_center + (i-1)*dx
    end

    # Fill faces
    start_face = start_center - dx/2
    @inbounds for i in 1:(total_zones + 1)
        faces[i] = start_face + (i-1)*dx
    end

    ilo = ghost_zones + 1
    ihi = ghost_zones + zones

    return UniformAxis{T}(
        centers,
        faces,
        dx,
        zones,
        ghost_zones,
        ilo,
        ihi,
        bounds
    )
end

#=
axis = UniformAxis((0.0, 2.0), 10, ghost_zones=3)

println(axis.dx)
println(axis.faces[axis.ilo-1:axis.ihi+1])
=#