




function apply_reflecting_boundaries!(U::ConservativeVariables, ng, nx)
    total = nx + 2ng 
    # Left boundary
    for i in 1:ng
        U.centers.density[ng - i + 1]        = U.centers.density[ng + i]
        U.centers.momentum[1][ng - i + 1]       = -U.centers.momentum[1][ng + i]
        U.centers.total_energy[ng - i + 1]   = U.centers.total_energy[ng + i]
    end

    # Right boundary
    for i in 1:ng
        U.centers.density[total - ng + i]       = U.centers.density[total - ng - i + 1]
        U.centers.momentum[1][total - ng + i]      = -U.centers.momentum[1][total - ng - i + 1]
        U.centers.total_energy[total - ng + i]  = U.centers.total_energy[total - ng - i + 1]
    end

    if U.faces !== nothing 
        total = length(U.faces.density)

        # Left boundary
        for i in 1:ng
            U.faces.density[ng - i + 1]        = U.faces.density[ng + i]
            U.faces.momentum[1][ng - i + 1]       = -U.faces.momentum[1][ng + i]
            U.faces.total_energy[ng - i + 1]   = U.faces.total_energy[ng + i]
        end

        # Right boundary
        for i in 1:ng
            U.faces.density[total - ng + i]       = U.faces.density[total - ng - i + 1]
            U.faces.momentum[1][total - ng + i]      = -U.faces.momentum[1][total - ng - i + 1]
            U.faces.total_energy[total - ng + i]  = U.faces.total_energy[total - ng - i + 1]
        end
    end

end

