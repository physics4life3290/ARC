



function Construct1DSodShockTubePrimitives(total_zones, ρL, uL, pL, ρR, uR, pR, γ)
    W = PrimitiveVariables(
        zeros(total_zones),
        zeros(total_zones),
        zeros(total_zones),
        zeros(total_zones)
    )
    mid = div(total_zones, 2)
    for i in 1:mid
        W.density[i] = ρL
        W.velocity[i] = uL
        W.pressure[i] = pL
        W.int_energy[i] = pL / (γ - 1) / ρL
    end
    for i in (mid+1):total_zones
        W.density[i] = ρR
        W.velocity[i] = uR
        W.pressure[i] = pR
        W.int_energy[i] = pR / (γ - 1) / ρR
    end
    return W
end