




function apply_periodic_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total = nx + 2ng

    # Fill left ghost zones (1:ng) with data from right physical zone (nx)
    for i in 1:ng
        U.density_centers[ng - i + 1]        = U.density_centers[nx + ng - i + 1]
        U.momentum_centers[ng - i + 1]       = U.momentum_centers[nx + ng - i + 1]
        U.total_energy_centers[ng - i + 1]   = U.total_energy_centers[nx + ng - i + 1]
    end

    # Fill right ghost zones (nx+ng+1 : total) with data from left physical zone (1:nx)
    for i in 1:ng
        U.density_centers[nx + ng + i]       = U.density_centers[ng + i]
        U.momentum_centers[nx + ng + i]      = U.momentum_centers[ng + i]
        U.total_energy_centers[nx + ng + i]  = U.total_energy_centers[ng + i]
    end
end