




function two_rarefaction_pressure(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64)

    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    numerator = cL + cR - 0.5 * (γ - 1) * (uR - uL)
    denominator = cL / pL^((γ - 1) / (2γ)) + cR / pR^((γ - 1) / (2γ))

    p_pv = (numerator / denominator)^(2γ / (γ - 1))

    return max(1e-8, p_pv) # positivity safeguard

end