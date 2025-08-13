




function Riemann_HLL(ρL, uL, pL, ρR, uR, pR, γ)
    # Toro's method - simplified pressure estimate
    p_star = max(1e-6, 0.5 * (pL + pR) - 0.125 * (uR - uL) * (ρL + ρR))
    u_star = 0.5 * (uL + uR) + 0.5 * (pL - pR) / (ρL + ρR)
    
    # Determine wave speeds (approximate)
    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)
    SL = minimum([uL - cL, u_star - sqrt(γ * p_star / ρL)])
    SR = maximum([uR + cR, u_star + sqrt(γ * p_star / ρR)])

    # Star region flux (HLL-type approach)
    if SL >= 0
        return flux(ρL, uL, pL)
    elseif SR <= 0
        return flux(ρR, uR, pR)
    else
        # HLL flux
        FL = flux(ρL, uL, pL)
        FR = flux(ρR, uR, pR)
        UL = conserved(ρL, uL, pL)
        UR = conserved(ρR, uR, pR)
        F1 = (SR * FL[1] - SL * FR[1] + SL * SR * (UR[1] - UL[1])) / (SR - SL)
        F2 = (SR * FL[2] - SL * FR[2] + SL * SR * (UR[2] - UL[2])) / (SR - SL)
        F3 = (SR * FL[3] - SL * FR[3] + SL * SR * (UR[3] - UL[3])) / (SR - SL)
        return F1, F2, F3
    end
end