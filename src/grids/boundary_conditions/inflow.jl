




function apply_inflow_boundaries!(U::ConservativeVariables, ng::Int, nx::Int;
                                  density_val::Float64=1.0,
                                  momentum_val::Float64=0.1,
                                  energy_val::Float64=1.0)
    total_centers = nx + 2ng

    # Left boundary inflow
    @inbounds @simd for i in 1:ng
        dest = ng - i + 1
        U.density_centers[dest]      = density_val
        U.momentum_centers[dest]     = momentum_val
        U.total_energy_centers[dest] = energy_val
    end

    # Right boundary inflow
    @inbounds @simd for i in 1:ng
        dest = total_centers - ng + i
        U.density_centers[dest]      = density_val
        U.momentum_centers[dest]     = momentum_val
        U.total_energy_centers[dest] = energy_val
    end

    return nothing
end
