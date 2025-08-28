




function apply_inflow_boundaries!(U::ConservativeVariables, ng::Int, nx::Int;
                                  density_val=1.0, momentum_val=0.1, energy_val=1.0)
    total = nx + 2ng

    # Left boundary inflow
    for i in 1:ng
        U.density_centers[ng - i + 1]      = density_val
        U.momentum_centers[ng - i + 1]     = momentum_val
        U.total_energy_centers[ng - i + 1] = energy_val
    end

    # Right boundary inflow (optional, could be zero or same)
    for i in 1:ng
        U.density_centers[total - ng + i]      = density_val
        U.momentum_centers[total - ng + i]     = momentum_val
        U.total_energy_centers[total - ng + i] = energy_val
    end
end
