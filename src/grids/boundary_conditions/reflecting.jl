



function apply_reflecting_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total_centers = nx + 2ng

    # Left boundary (centers)
    @inbounds @simd for i in 1:ng
        U.density_centers[i]      = U.density_centers[ng + i]
        U.momentum_centers[i]     = -U.momentum_centers[ng + i]
        U.total_energy_centers[i] = U.total_energy_centers[ng + i]
    end

    # Right boundary (centers)
    @inbounds @simd for i in 1:ng
        src = nx + ng - i + 1
        dest = nx + ng + i
        U.density_centers[dest]      = U.density_centers[src]
        U.momentum_centers[dest]     = -U.momentum_centers[src]
        U.total_energy_centers[dest] = U.total_energy_centers[src]
    end

    # Faces (if they exist)
    if U.density_faces !== nothing
        total_faces = length(U.density_faces)

        # Left boundary (faces)
        @inbounds @simd for i in 1:ng
            U.density_faces[i]      = U.density_faces[ng + i]
            U.momentum_faces[i]     = -U.momentum_faces[ng + i]
            U.total_energy_faces[i] = U.total_energy_faces[ng + i]
        end

        # Right boundary (faces)
        @inbounds @simd for i in 1:ng
            src = total_faces - ng - i + 1
            dest = total_faces - ng + i
            U.density_faces[dest]      = U.density_faces[src]
            U.momentum_faces[dest]     = -U.momentum_faces[src]
            U.total_energy_faces[dest] = U.total_energy_faces[src]
        end
    end

    return nothing
end

