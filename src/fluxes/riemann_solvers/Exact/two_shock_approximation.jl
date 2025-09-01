




function two_shock_pressure(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64)

    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    # Shock impedance factors (Eq. 4.47)
    gL = sqrt((2 / ((γ + 1) * ρL)) / (pL + (γ - 1) / (γ + 1) * pL))
    gR = sqrt((2 / ((γ + 1) * ρR)) / (pR + (γ - 1) / (γ + 1) * pR))

    numerator = gL * pL + gR * pR - (uR - uL)
    denominator = gL + gR

    p_ts = numerator / denominator

    return max(1e-8, p_ts) # positivity safeguard
    
end