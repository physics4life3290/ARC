




function apply_periodic_boundaries!(U::ConservativeVariables, total::Int32)
    

    # Fill left ghost zones (1:ng) with data from right physical zone (nx)
    for i in 1:ng
        U.centers.density[ng - i + 1]        = U.centers.density[nx + ng - i + 1]
        U.centers.momentum[1][ng - i + 1]       = U.centers.momentum[1][nx + ng - i + 1]
        U.centers.total_energy[ng - i + 1]   = U.centers.total_energy[nx + ng - i + 1]
    end

    # Fill right ghost zones (nx+ng+1 : total) with data from left physical zone (1:nx)
    for i in 1:ng
        U.centers.density[nx + ng + i]       = U.centers.density[ng + i]
        U.centers.momentum[1][nx + ng + i]      = U.centers.momentum[1][ng + i]
        U.centers.total_energy[nx + ng + i]  = U.centers.total_energy[ng + i]
    end
end