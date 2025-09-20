




function apply_robin_boundaries!(U::ConservativeVariables, ng::Int, nx::Int;
                                 alpha::Float64=1.0, beta::Float64=1.0,
                                 gamma_density::Float64=1.0,
                                 gamma_momentum::Float64=0.0,
                                 gamma_energy::Float64=1.0)

    total = nx + 2ng

    # Left boundary
    @inbounds @simd for i in 1:ng
        inside = ng + i
        dest = ng - i + 1

        U.density_centers[dest]      = (gamma_density - alpha * U.density_centers[inside]) / beta
        U.momentum_centers[dest]     = (gamma_momentum - alpha * U.momentum_centers[inside]) / beta
        U.total_energy_centers[dest] = (gamma_energy - alpha * U.total_energy_centers[inside]) / beta
    end

    # Right boundary
    @inbounds @simd for i in 1:ng
        inside = total - ng - i + 1
        dest   = total - ng + i

        U.density_centers[dest]      = (gamma_density - alpha * U.density_centers[inside]) / beta
        U.momentum_centers[dest]     = (gamma_momentum - alpha * U.momentum_centers[inside]) / beta
        U.total_energy_centers[dest] = (gamma_energy - alpha * U.total_energy_centers[inside]) / beta
    end

    return nothing
end
