




function apply_outflow_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total = nx + 2ng

    # Left boundary: copy innermost left values into ghost zones
    for i in 1:ng
        U.density_centers[ng - i + 1]        = U.density_centers[ng + 1]
        U.momentum_centers[ng - i + 1]       = U.momentum_centers[ng + 1]
        U.total_energy_centers[ng - i + 1]   = U.total_energy_centers[ng + 1]
    end

    # Right boundary: copy innermost right values into ghost zones
    for i in 1:ng
        U.density_centers[nx + ng + i]       = U.density_centers[nx + ng]
        U.momentum_centers[nx + ng + i]      = U.momentum_centers[nx + ng]
        U.total_energy_centers[nx + ng + i]  = U.total_energy_centers[nx + ng]
    end
end