




struct LogarithmicAxis{T<:AbstractFloat}
    centers::Vector{T}
    faces::Vector{T}
    dx::Vector{T}
    zones::Int
    ghost_zones::Int
    bounds::Tuple{T,T}
end

function SymLogAxis(domain::Tuple{T}, zones::Int; ghost_zones::Int=3) where T
    @assert zones > 0 "zones must be positive"
    @assert ghost_zones ≥ 0 "ghost_zones must be non-negative"

    xmax = domain[1]/2
    xmin = -xmax
    bounds = (xmin, xmax)

    # integer division for positive/negative halves
    n_half = div(zones, 2)

    # positive interior centers (log-spaced)
    pos_centers = 10 .^ range(log10(1e-8), log10(xmax), length=n_half)
    neg_centers = -reverse(pos_centers)
    interior_centers = vcat(neg_centers, pos_centers)

    # compute dx of first and last interior cells
    dx_left  = interior_centers[2] - interior_centers[1]
    dx_right = interior_centers[end] - interior_centers[end-1]

    # prepend ghost centers on left, append on right
    ghost_left  = interior_centers[1] .- dx_left .* (ghost_zones:-1:1)
    ghost_right = interior_centers[end] .+ dx_right .* (1:ghost_zones)
    centers = vcat(ghost_left, interior_centers, ghost_right)

    # faces are edges between centers
    faces = vcat(
        centers[1] - dx_left/2,           # leftmost face
        (centers[1:end-1] .+ centers[2:end]) ./ 2,  # midpoints
        centers[end] + dx_right/2         # rightmost face
    )

    dx = diff(faces)  # cell widths

    return LogarithmicAxis(centers, faces, dx, zones, ghost_zones, bounds)
end




