




function apply_reflecting_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total = nx + 2ng

    # Left boundary
    for i in 1:ng
        U.density_centers[ng - i + 1]        = U.density_centers[ng + i]
        U.momentum_centers[ng - i + 1]       = -U.momentum_centers[ng + i]
        U.total_energy_centers[ng - i + 1]   = U.total_energy_centers[ng + i]
    end

    # Right boundary
    for i in 1:ng
        U.density_centers[total - ng + i]       = U.density_centers[total - ng - i + 1]
        U.momentum_centers[total - ng + i]      = -U.momentum_centers[total - ng - i + 1]
        U.total_energy_centers[total - ng + i]  = U.total_energy_centers[total - ng - i + 1]
    end
end

