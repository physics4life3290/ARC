




function apply_robin_boundaries!(U::ConservativeVariables, ng::Int, nx::Int;
                                 alpha=1.0, beta=1.0, gamma_density=1.0,
                                 gamma_momentum=0.0, gamma_energy=1.0)
    total = nx + 2ng

    # Helper function
    robin!(val_inside, val_out, γ) = (γ - alpha * val_inside) / beta

    # Left boundary
    for i in 1:ng
        U.density_centers[ng - i + 1]      = robin!(U.density_centers[ng + i], U.density_centers[ng - i + 1], gamma_density)
        U.momentum_centers[ng - i + 1]     = robin!(U.momentum_centers[ng + i], U.momentum_centers[ng - i + 1], gamma_momentum)
        U.total_energy_centers[ng - i + 1] = robin!(U.total_energy_centers[ng + i], U.total_energy_centers[ng - i + 1], gamma_energy)
    end

    # Right boundary
    for i in 1:ng
        U.density_centers[total - ng + i]      = robin!(U.density_centers[total - ng - i + 1], U.density_centers[total - ng + i], gamma_density)
        U.momentum_centers[total - ng + i]     = robin!(U.momentum_centers[total - ng - i + 1], U.momentum_centers[total - ng + i], gamma_momentum)
        U.total_energy_centers[total - ng + i] = robin!(U.total_energy_centers[total - ng - i + 1], U.total_energy_centers[total - ng + i], gamma_energy)
    end
end
