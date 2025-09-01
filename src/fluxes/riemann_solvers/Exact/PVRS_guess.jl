




function PVRS_guess(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64)

    aL, aR = cs(γ, pL, ρL), cs(γ, pR, ρR)
    p̃ = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (ρL + ρR) * (aL + aR)
    p0 = max(1e-8, p̃)

    return p0

end