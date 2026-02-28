




include("grid_structs.jl")

# ----------------------
# Factory function
# ----------------------
function make_grid(kind::Symbol, origin::NTuple{N,T}, domain::NTuple{N,T},
                   nx::NTuple{N,Int32}, ghost_zones::Int32;
                   threshold::Float64=0.02, dx_array=nothing) where {T,N}

    if kind === :uniform
        axes = ntuple(i -> begin
            dx = (domain[i] - origin[i]) / nx[i]

            # include ghost zones
            centers = collect(origin[i] .+ dx*(0.5:(nx[i]-0.5)))
            centers = vcat(centers[1] .- dx*reverse(1:ghost_zones), centers, centers[end] .+ dx*(1:ghost_zones))

            faces = collect(origin[i] .+ dx*(0:nx[i]))
            faces = vcat(faces[1] .- dx*reverse(1:ghost_zones), faces, faces[end] .+ dx*(1:ghost_zones))

            UniformAxis(centers, faces, dx, Int32(length(centers)), (origin[i], origin[i]+domain[i]))
        end, N)

        active_zones = ntuple(i -> (ghost_zones+1):(ghost_zones+nx[i]), N)
        total_zones = ntuple(i -> Int32(nx[i] + 2*ghost_zones), N)

        return UniformGrid{T,N}(axes, origin, ghost_zones, total_zones, active_zones, domain)

    elseif kind === :adaptive
        root_blocks = [Block{T}(Vector{Vector{T}}(), zero(T), 0, Block{T}[], origin[1], nx[1])]
        return AdaptiveGrid{T,N}(root_blocks, origin, domain, ghost_zones, 0)

    else
        error("Unknown grid type: $kind")
    end
end

abstract type Verbosity end
struct Standard <: Verbosity end
struct Verbose  <: Verbosity end

log(::Standard, args...) = nothing
log(::Verbose, args...)  = println(args...)

function build_uniform_axis(i, origin, domain, nx, ghost_zones, ::Type{T}, mode::Verbosity) where T

    log(mode, "\n→ AXIS ", i, " INITIALIZATION")

    dx = (domain[i] - origin[i]) / nx[i]
    log(mode, "   Δx = ", dx)

    centers = origin[i] .+ dx*(0.5:(nx[i]-0.5))
    centers = vcat(
        centers[1] .- dx*reverse(1:ghost_zones),
        centers,
        centers[end] .+ dx*(1:ghost_zones)
    )

    faces = origin[i] .+ dx*(0:nx[i])
    faces = vcat(
        faces[1] .- dx*reverse(1:ghost_zones),
        faces,
        faces[end] .+ dx*(1:ghost_zones)
    )

    log(mode, "   Final center count: ", length(centers))
    log(mode, "   Final face count:   ", length(faces))

    return UniformAxis(
        collect(centers),
        collect(faces),
        dx,
        Int32(length(centers)),
        (origin[i], origin[i] + domain[i])
    )
end

function make_grid(kind::Symbol, origin::NTuple{N,T}, domain::NTuple{N,T},
                   nx::NTuple{N,Int32}, ghost_zones::Int32;
                   threshold=0.02, dx_array=nothing,
                   mode::Verbosity=Standard()) where {T,N}

    log(mode, "══════════════════════════════════════════════")
    log(mode, "GRID CONSTRUCTION START")

    if kind === :uniform

        axes = ntuple(i -> build_uniform_axis(i, origin, domain, nx, ghost_zones, T, mode), N)

        active_zones = ntuple(i -> (ghost_zones+1):(ghost_zones+nx[i]), N)
        total_zones  = ntuple(i -> Int32(nx[i] + 2*ghost_zones), N)

        return UniformGrid{T,N}(axes, origin, ghost_zones, total_zones, active_zones, domain)

    elseif kind === :adaptive

        root_blocks = [Block{T}(Vector{Vector{T}}(), zero(T), 0, Block{T}[], origin[1], nx[1])]
        return AdaptiveGrid{T,N}(root_blocks, origin, domain, ghost_zones, 0)

    else
        error("Unknown grid type: $kind")
    end
end