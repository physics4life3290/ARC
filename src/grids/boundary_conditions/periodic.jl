




function apply_periodic_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total_centers = nx + 2ng

    # Left ghost zones (1:ng) ← data from right physical zone
    @inbounds @simd for i in 1:ng
        src = nx + ng - i + 1
        dest = ng - i + 1
        U.density_centers[dest]      = U.density_centers[src]
        U.momentum_centers[dest]     = U.momentum_centers[src]
        U.total_energy_centers[dest] = U.total_energy_centers[src]
    end

    # Right ghost zones (nx+ng+1 : total) ← data from left physical zone
    @inbounds @simd for i in 1:ng
        src = ng + i
        dest = nx + ng + i
        U.density_centers[dest]      = U.density_centers[src]
        U.momentum_centers[dest]     = U.momentum_centers[src]
        U.total_energy_centers[dest] = U.total_energy_centers[src]
    end

    return nothing
end
