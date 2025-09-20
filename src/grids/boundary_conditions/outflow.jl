




function apply_outflow_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    # Left ghost zones: copy innermost left value
    left_val_idx = ng + 1
    @inbounds @simd for i in 1:ng
        dest = ng - i + 1
        U.density_centers[dest]      = U.density_centers[left_val_idx]
        U.momentum_centers[dest]     = U.momentum_centers[left_val_idx]
        U.total_energy_centers[dest] = U.total_energy_centers[left_val_idx]
    end

    # Right ghost zones: copy innermost right value
    right_val_idx = nx + ng
    @inbounds @simd for i in 1:ng
        dest = nx + ng + i
        U.density_centers[dest]      = U.density_centers[right_val_idx]
        U.momentum_centers[dest]     = U.momentum_centers[right_val_idx]
        U.total_energy_centers[dest] = U.total_energy_centers[right_val_idx]
    end

    return nothing
end
