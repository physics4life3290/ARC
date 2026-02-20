




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

    elseif kind === :nonuniform
        if dx_array === nothing
            error("dx_array must be supplied for nonuniform grids")
        end

        axes = ntuple(i -> begin
            dx_i = dx_array[i]

            centers = origin[i] .+ cumsum(dx_i) .- 0.5*dx_i
            faces = [origin[i]; origin[i] .+ cumsum(dx_i)]

            # ghost zones: repeat first/last dx
            pre_centers = centers[1] .- reverse(cumsum(dx_i[1:ghost_zones]))
            post_centers = centers[end] .+ cumsum(dx_i[end-ghost_zones+1:end])
            centers_ext = vcat(pre_centers, centers, post_centers)

            pre_faces = faces[1] .- reverse(cumsum(dx_i[1:ghost_zones]))
            post_faces = faces[end] .+ cumsum(dx_i[end-ghost_zones+1:end])
            faces_ext = vcat(pre_faces, faces, post_faces)

            NonUniformAxis(centers_ext, faces_ext, dx_i, Int32(length(centers_ext)), (origin[i], origin[i]+domain[i]))
        end, N)

        active_zones = ntuple(i -> (ghost_zones+1):(ghost_zones+nx[i]), N)
        total_zones = ntuple(i -> Int32(nx[i] + 2*ghost_zones), N)

        return NonUniformGrid{T,N}(axes, origin, ghost_zones, total_zones, active_zones, domain)

    elseif kind === :adaptive
        root_blocks = [Block{T}(Vector{Vector{T}}(), zero(T), 0, Block{T}[], origin[1], nx[1])]
        return AdaptiveGrid{T,N}(root_blocks, origin, domain, ghost_zones, 0)

    else
        error("Unknown grid type: $kind")
    end
end

