




function apply_outflow_boundaries!(U::ConservativeVariables, ng::Int32, nx::Int32)
    total = nx + 2ng

    # Left boundary: copy innermost left values into ghost zones
    for i in 1:ng
        U.centers.density[ng - i + 1]        = U.centers.density[ng + 1]
        U.centers.momentum[1][ng - i + 1]       = U.centers.momentum[1][ng + 1]
        U.centers.total_energy[ng - i + 1]   = U.centers.total_energy[ng + 1]
    end

    # Right boundary: copy innermost right values into ghost zones
    for i in 1:ng
        U.centers.density[nx + ng + i]       = U.centers.density[nx + ng]
        U.centers.momentum[1][nx + ng + i]      = U.centers.momentum[1][nx + ng]
        U.centers.total_energy[nx + ng + i]  = U.centers.total_energy[nx + ng]
    end
end