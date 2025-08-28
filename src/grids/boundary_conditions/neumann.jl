




function apply_neumann_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total = nx + 2ng

    # Left boundary
    for i in 1:ng
        U.density_centers[ng - i + 1]      = U.density_centers[ng + 1]
        U.momentum_centers[ng - i + 1]     = U.momentum_centers[ng + 1]
        U.total_energy_centers[ng - i + 1] = U.total_energy_centers[ng + 1]
    end

    # Right boundary
    for i in 1:ng
        U.density_centers[total - ng + i]      = U.density_centers[total - ng]
        U.momentum_centers[total - ng + i]     = U.momentum_centers[total - ng]
        U.total_energy_centers[total - ng + i] = U.total_energy_centers[total - ng]
    end

    if U.density_faces !== nothing
        total = length(U.density_faces)

        for i in 1:ng
            U.density_faces[ng - i + 1]      = U.density_faces[ng + 1]
            U.momentum_faces[ng - i + 1]     = U.momentum_faces[ng + 1]
            U.total_energy_faces[ng - i + 1] = U.total_energy_faces[ng + 1]
        end
        for i in 1:ng
            U.density_faces[total - ng + i]      = U.density_faces[total - ng]
            U.momentum_faces[total - ng + i]     = U.momentum_faces[total - ng]
            U.total_energy_faces[total - ng + i] = U.total_energy_faces[total - ng]
        end
    end
end
